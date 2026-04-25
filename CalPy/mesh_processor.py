import gmsh
import re

# ── SETTINGS ────────────────────────────────────────────────────────────────
input_msh        = "UCRMWS.msh"       # Pointwise .msh input file
output_inp       = "mesh_C3D4.nam"       # Final cleaned Abaqus .inp
skin_elset_inp   = "skinElementList_c3D4.nam"    # Output ELSET file for skin elements
skin_elset_name  = "SkinElement"            # Name of the skin ELSET
skin_nset        = "Surface"                # NSET to define skin elements against

# Set to True to upgrade C3D4 -> C3D10, False to keep as C3D4
upgrade_to_C3D10 = False

# NSETs to drop from the final output
drop_nsets = {'EALL', 'UNSPECIFIED'}
# ────────────────────────────────────────────────────────────────────────────


def upgrade_mesh(input_msh, temp_inp, upgrade):
    """
    Load .msh from Pointwise, optionally upgrade to second order,
    export raw .inp via Gmsh API.
    """
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)

    gmsh.merge(input_msh)

    print("\nElement types BEFORE upgrade:")
    for dim, tag in gmsh.model.getEntities():
        for et in gmsh.model.mesh.getElements(dim, tag)[0]:
            info = gmsh.model.mesh.getElementProperties(et)
            print(f"  dim={dim} tag={tag} type={et} ({info[0]}) nodesPerElem={info[3]}")

    if upgrade:
        # Note: on an imported mesh with no CAD geometry, Gmsh places midside
        # nodes at edge midpoints (equivalent to SecondOrderLinear = 1)
        gmsh.model.mesh.setOrder(2)
        print("\nElement types AFTER upgrade:")
        for dim, tag in gmsh.model.getEntities():
            for et in gmsh.model.mesh.getElements(dim, tag)[0]:
                info = gmsh.model.mesh.getElementProperties(et)
                print(f"  dim={dim} tag={tag} type={et} ({info[0]}) nodesPerElem={info[3]}")
    else:
        print("\nSkipping upgrade — keeping linear elements (C3D4)")

    gmsh.option.setNumber("Mesh.SaveGroupsOfNodes",    2)
    gmsh.option.setNumber("Mesh.SaveGroupsOfElements", 2)
    gmsh.option.setNumber("Mesh.SaveAll", 0)

    gmsh.write(temp_inp)
    gmsh.finalize()
    print(f"\nRaw Gmsh export written to: {temp_inp}")


def clean_inp(input_file, output_file, drop_nsets, upgrade):
    """
    Post-process raw Gmsh .inp:
      - Keep *NODE block
      - Keep *ELEMENT block (C3D10 or C3D4), rename ELSET to EALL
      - Drop all surface ELEMENTs (CPS6 or CPS3)
      - Keep *NSET blocks except those in drop_nsets
    Returns the parsed elements and nsets for downstream use.
    """
    with open(input_file, 'r') as f:
        blocks = re.split(r'(?=^\*)', f.read(), flags=re.MULTILINE)

    node_block        = None
    volume_elem_block = None
    keep_nsets        = {}
    elements          = {}  # elem_id -> [node ids]
    nsets             = {}  # nset_name -> set of node ids

    target_elem_type = 'C3D10' if upgrade else 'C3D4'

    for block in blocks:
        lines = block.strip().splitlines()
        if not lines:
            continue
        header    = lines[0]
        header_up = header.upper()

        # NODE block
        if re.match(r'^\*NODE', header_up):
            header_line = lines[0]
            rounded = [header_line]
            for line in lines[1:]:
                parts = line.split(',')
                if len(parts) >= 4:
                    node_id = parts[0].strip()
                    coords  = [f'{float(p):.7f}' for p in parts[1:4]]
                    rounded.append(f'  {node_id}, ' + ', '.join(coords))
                else:
                    rounded.append(line)
            node_block = rounded

        # Volume element block — C3D10 or C3D4 depending on upgrade flag
        elif re.match(r'^\*ELEMENT', header_up) and target_elem_type in header_up:
            header_line = f'*ELEMENT, type={target_elem_type}, ELSET=EALL'
            reindexed = [header_line]
            new_idx = 1
            for line in lines[1:]:
                parts = [p.strip() for p in line.split(',')]
                if parts and parts[0].isdigit():
                    nodes = ', '.join(parts[1:])
                    reindexed.append(f'{new_idx}, {nodes}')
                    # Parse with new index for skin extraction
                    node_ids = [int(n) for n in parts[1:] if n]
                    elements[new_idx] = node_ids
                    new_idx += 1
            volume_elem_block = reindexed

        # Surface elements (CPS6, CPS3) — drop
        elif re.match(r'^\*ELEMENT', header_up):
            pass

        # NSET blocks
        elif re.match(r'^\*NSET', header_up):
            match = re.search(r'NSET\s*=\s*(\w+)', header, re.IGNORECASE)
            if match:
                name     = match.group(1)
                node_ids = set()
                for line in lines[1:]:
                    for token in line.split(','):
                        token = token.strip()
                        if token.isdigit():
                            node_ids.add(int(token))
                nsets[name] = node_ids
                if name.upper() not in drop_nsets:
                    keep_nsets[name] = lines

    if node_block is None:
        raise ValueError("No *NODE block found")
    if volume_elem_block is None:
        raise ValueError(
            f"No {target_elem_type} *ELEMENT block found — "
            f"{'upgrade may have failed' if upgrade else 'check your mesh'}"
        )

    with open(output_file, 'w') as f:
        f.write('*Heading\n')
        f.write(f' {output_file}\n\n')
        f.write('\n'.join(node_block) + '\n\n')
        f.write('\n'.join(volume_elem_block) + '\n\n')
        for name, lines in keep_nsets.items():
            f.write('\n'.join(lines) + '\n\n')

    print(f"\nCleaned output written to: {output_file}")
    print(f"  ELSET : EALL ({target_elem_type})")
    print(f"  NSETs : {list(keep_nsets.keys())}")

    return elements, nsets


def write_skin_elset(elements, nsets, target_nset, elset_name, output_file):
    """
    Find all elements that share at least one node with target_nset,
    write them as an ELSET to output_file.
    """
    if target_nset not in nsets:
        raise ValueError(
            f"NSET '{target_nset}' not found. "
            f"Available NSETs: {list(nsets.keys())}"
        )

    nset_nodes  = nsets[target_nset]
    skin_elems  = sorted(
        elem_id for elem_id, nodes in elements.items()
        if any(n in nset_nodes for n in nodes)
    )

    with open(output_file, 'w') as f:
        #f.write(f'*ELSET, ELSET={elset_name}\n')
        for i in range(0, len(skin_elems)):
            f.write(str(skin_elems[i]) + '\n')

    print(f"\nSkin ELSET written to: {output_file}")
    print(f"  ELSET        : {elset_name}")
    print(f"  Source NSET  : {target_nset} ({len(nset_nodes)} nodes)")
    print(f"  Elements     : {len(skin_elems)}")


if __name__ == "__main__":
    temp_inp = "mesh_raw.inp"

    # Step 1 — Load, optionally upgrade, export raw .inp
    upgrade_mesh(input_msh, temp_inp, upgrade_to_C3D10)

    # Step 2 — Clean raw .inp, return parsed elements and nsets
    elements, nsets = clean_inp(temp_inp, output_inp, drop_nsets, upgrade_to_C3D10)

    # Step 3 — Write skin ELSET based on target NSET
    write_skin_elset(elements, nsets, skin_nset, skin_elset_name, skin_elset_inp)
