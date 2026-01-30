#!/usr/bin/env python3
import argparse
import os


def read_su2_markers(su2_file):
    """
    Reads an SU2 mesh file and returns:
    { marker_name : sorted list of unique node indices (1-based) }
    """
    markers = {}
    current_marker = None
    elems_remaining = 0

    with open(su2_file, "r") as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("%"):
                continue

            if line.startswith("MARKER_TAG"):
                current_marker = line.split("=")[1].strip()
                markers[current_marker] = set()

            elif line.startswith("MARKER_ELEMS"):
                elems_remaining = int(line.split("=")[1].strip())

            elif elems_remaining > 0 and current_marker is not None:
                parts = line.split()
                # parts[0] = element type
                node_ids = [int(p) + 1 for p in parts[1:]]  # 0-based → 1-based
                markers[current_marker].update(node_ids)
                elems_remaining -= 1

    return {k: sorted(v) for k, v in markers.items()}


def read_su2_elements(su2_file):
    """
    Reads element connectivity from SU2 mesh.
    Returns list of elements where each element is a dict with:
    { 'type': elem_type, 'nodes': [node_ids in 1-based], 'id': element_id (1-based) }
    """
    elements = []
    in_elements = False

    with open(su2_file, "r") as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("%"):
                continue

            if line.startswith("NELEM"):
                in_elements = True
                continue

            if in_elements:
                if line.startswith("NPOIN") or line.startswith("NMARK"):
                    break

                parts = line.split()
                elem_type = int(parts[0])
                # For C3D4 (type 10), we have 4 nodes + element ID
                node_ids = [int(p) + 1 for p in parts[1:-1]]  # Convert to 1-based
                elem_id = int(parts[-1]) + 1  # Convert to 1-based

                elements.append({
                    'type': elem_type,
                    'nodes': node_ids,
                    'id': elem_id
                })

    return elements


def read_su2_pressure_marker(su2_file, marker_name):
    """
    Reads boundary elements for a specific pressure marker.
    Returns list of boundary faces, where each face is a list of node IDs (1-based).
    """
    faces = []
    current_marker = None
    elems_remaining = 0

    with open(su2_file, "r") as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("%"):
                continue

            if line.startswith("MARKER_TAG"):
                current_marker = line.split("=")[1].strip()

            elif line.startswith("MARKER_ELEMS"):
                elems_remaining = int(line.split("=")[1].strip())

            elif elems_remaining > 0 and current_marker == marker_name:
                parts = line.split()
                # parts[0] = element type (5 for triangle)
                node_ids = [int(p) + 1 for p in parts[1:]]  # Convert to 1-based
                faces.append(node_ids)
                elems_remaining -= 1

    return faces


def find_opposite_node_c3d4(tet_nodes, face_nodes):
    """
    For a C3D4 tetrahedron, given 3 face nodes, find the opposite node.
    Returns the P-number (1-4) based on CalculiX convention.
    
    CalculiX C3D4 face definitions for element with nodes [n1, n2, n3, n4]:
    - P1 (S1): face [n1, n2, n3] - opposite to n4 (position 4, index 3)
    - P2 (S2): face [n1, n4, n2] - opposite to n3 (position 3, index 2)
    - P3 (S3): face [n2, n4, n3] - opposite to n1 (position 1, index 0)
    - P4 (S4): face [n3, n4, n1] - opposite to n2 (position 2, index 1)
    """
    face_set = set(face_nodes)
    tet_set = set(tet_nodes)

    # Check if face nodes are subset of tet nodes
    if not face_set.issubset(tet_set):
        return None

    # Find the node in tet that's not in face
    opposite_nodes = tet_set - face_set
    if len(opposite_nodes) != 1:
        return None

    opposite_node = list(opposite_nodes)[0]

    # Map opposite node index to P number based on CalculiX C3D4 convention
    try:
        opposite_idx = tet_nodes.index(opposite_node)
        if opposite_idx == 3:
            return 1  # P1: opposite to 4th node
        elif opposite_idx == 2:
            return 2  # P2: opposite to 3rd node
        elif opposite_idx == 0:
            return 3  # P3: opposite to 1st node
        elif opposite_idx == 1:
            return 4  # P4: opposite to 2nd node
    except ValueError:
        return None

    return None


def process_pressure_bc(su2_file, marker_name):
    """
    Process pressure BC marker and return element-to-P mapping.
    Returns dict: { P_number: [list of element IDs (1-based)] }
    
    Optimized version using hash map for O(n+m) complexity instead of O(n*m).
    """
    elements = read_su2_elements(su2_file)
    faces = read_su2_pressure_marker(su2_file, marker_name)

    # Build a hash map: frozenset(3 nodes) -> list of (element_id, p_number)
    # This allows O(1) lookup instead of iterating through all elements for each face
    face_to_elem = {}
    
    for elem in elements:
        if elem['type'] != 10:  # Only C3D4
            continue
        
        # Generate all 4 possible faces for this tetrahedral element
        # and store which P-number each face corresponds to
        nodes = elem['nodes']
        
        # P1: opposite to 4th node (index 3) -> face is nodes [0,1,2]
        face_p1 = frozenset([nodes[0], nodes[1], nodes[2]])
        # P2: opposite to 3rd node (index 2) -> face is nodes [0,3,1]
        face_p2 = frozenset([nodes[0], nodes[3], nodes[1]])
        # P3: opposite to 1st node (index 0) -> face is nodes [1,3,2]
        face_p3 = frozenset([nodes[1], nodes[3], nodes[2]])
        # P4: opposite to 2nd node (index 1) -> face is nodes [2,3,0]
        face_p4 = frozenset([nodes[2], nodes[3], nodes[0]])
        
        face_to_elem[face_p1] = (elem['id'], 1)
        face_to_elem[face_p2] = (elem['id'], 2)
        face_to_elem[face_p3] = (elem['id'], 3)
        face_to_elem[face_p4] = (elem['id'], 4)

    # Now process boundary faces with O(1) lookup
    p_mapping = {1: [], 2: [], 3: [], 4: []}
    
    for face in faces:
        face_key = frozenset(face)
        if face_key in face_to_elem:
            elem_id, p_num = face_to_elem[face_key]
            p_mapping[p_num].append(elem_id)

    return p_mapping


def write_nam_files(markers, output_dir, skip_markers=None):
    """
    Writes one .nam file per marker.
    Skips marker named 'Unspecified' and any markers in skip_markers list.
    
    Args:
        markers: Dictionary of marker names to node sets
        output_dir: Directory to write .nam files
        skip_markers: List of marker names to skip (e.g., pressure BC markers)
    """
    os.makedirs(output_dir, exist_ok=True)
    
    if skip_markers is None:
        skip_markers = []

    for marker, nodes in markers.items():
        if marker.lower() == "unspecified":
            continue
        
        if marker in skip_markers:
            continue

        filepath = os.path.join(output_dir, f"{marker}.nam")
        with open(filepath, "w") as f:
            f.write(f"*NSET, NSET={marker}\n")
            for n in nodes:
                f.write(f"{n}\n")


def write_pressure_nam_files(p_mapping, marker_name, output_dir):
    """
    Writes .nam files for pressure BC with element IDs and P numbers.
    Only writes files for P numbers that have elements.
    """
    os.makedirs(output_dir, exist_ok=True)

    for p_num in [1, 2, 3, 4]:
        if not p_mapping[p_num]:
            continue  # Skip if no elements for this P number

        filepath = os.path.join(output_dir, f"{marker_name}P{p_num}.nam")
        with open(filepath, "w") as f:
            f.write(f"*ELSET, ELSET={marker_name}P{p_num}\n")
            for elem_id in p_mapping[p_num]:
                f.write(f"{elem_id},\n")


def main():
    parser = argparse.ArgumentParser(
        description="Extract marker node sets from an SU2 mesh and write .nam files."
    )
    parser.add_argument(
        "--mesh",
        required=True,
        help="Path to the SU2 mesh file (e.g. mesh.su2)"
    )
    parser.add_argument(
        "--out",
        default="markers",
        help="Output directory for .nam files (default: ./markers)"
    )
    parser.add_argument(
        "--pBC",
        help="Name of pressure BC marker to process for element-based DLOAD"
    )

    args = parser.parse_args()

    if not os.path.isfile(args.mesh):
        raise FileNotFoundError(f"SU2 mesh file not found: {args.mesh}")

    # Always process regular markers, but skip the pBC marker if specified
    markers = read_su2_markers(args.mesh)
    skip_markers = [args.pBC] if args.pBC else []
    write_nam_files(markers, args.out, skip_markers)

    num_written = len([m for m in markers if m.lower() != 'unspecified' and m not in skip_markers])
    print(f"Processed {args.mesh}")
    print(f"Wrote {num_written} marker sets to '{args.out}/'")

    # Process pressure BC if specified
    if args.pBC:
        p_mapping = process_pressure_bc(args.mesh, args.pBC)
        write_pressure_nam_files(p_mapping, args.pBC, args.out)

        num_p_files = sum(1 for v in p_mapping.values() if v)
        print(f"Wrote {num_p_files} pressure BC files for marker '{args.pBC}'")
        for p_num in [1, 2, 3, 4]:
            if p_mapping[p_num]:
                print(f"  {args.pBC}P{p_num}.nam: {len(p_mapping[p_num])} elements")


if __name__ == "__main__":
    main()