#!/usr/bin/env python3

"""
Generate an initial density field for topology optimization while enforcing
passive (solid) elements.

This script creates a file called `density.dat` containing one density value
per design variable (typically one per finite element). All elements are
initialized to a prescribed default density value (e.g., 0.3). Elements listed
in a provided `.nam` file are treated as passive elements and are forced to
have density equal to 1.0.

The `.nam` file is expected to contain one element ID per line. Element IDs
are assumed to be **1-based**, consistent with CalculiX element numbering.
Internally, the script converts these IDs to Python's **0-based indexing**
before modifying the design variable array.

Logic follows the same behavior as the C function:

    design[eid - 1] = 1.0

where `eid` is a valid element ID satisfying:

    1 <= eid <= number_of_design_variables

If an element ID falls outside this range, a warning is printed.

Outputs
-------
density.dat
    Text file containing one density value per line with 15 decimal places.
    This file can be used as the initial density field for CalTop or other
    topology optimization workflows.

Inputs
------
skinElementList.nam
    File containing 1-based element IDs that must remain solid (passive
    elements). One ID per line.

Parameters
----------
n : int
    Number of design variables (typically equal to the number of elements).
rho_default : float
    Initial density assigned to all elements before passive elements are
    enforced.

Workflow
--------
1. Initialize the density array with the prescribed default density.
2. Read passive element IDs from the `.nam` file.
3. For each valid passive element ID:
       design[eid - 1] = 1.0
4. Write the resulting density field to `density.dat`.

This script is intended for preprocessing topology optimization problems
using the CalTop framework.
"""

# Number of design variables
n = 1485336

# Default density for all elements
rho_default = 0.1

# Input file containing 1-based element IDs, one per line
passive_file = "skinElementList.nam"

# Output density file
output_file = "density.dat"


def read_passive_ids(filename):
    passive_ids = []

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()

            # Skip blank lines
            if not line:
                continue

            # Read integer element ID
            eid = int(line)
            passive_ids.append(eid)

    return passive_ids


def filter_out_passive_elems_density(design, ne, passive_ids):
    for eid in passive_ids:   # eid is 1-based
        if 1 <= eid <= ne:
            design[eid - 1] = 1.0   # convert to 0-based index
        else:
            print(f"Warning: Passive element ID {eid} is out of bounds [1, {ne}]")


def write_density_file(filename, design):
    with open(filename, "w") as f:
        for rho in design:
            f.write(f"{rho:.15f}\n")


def main():
    # Initialize all design variables to the prescribed density
    design = [rho_default] * n

    # Read passive element IDs from file
    passive_ids = read_passive_ids(passive_file)

    # Apply same logic as the C function
    filter_out_passive_elems_density(design, n, passive_ids)

    # Write final density file
    write_density_file(output_file, design)

    print(f"{output_file} written with {n} entries.")
    print(f"Default density = {rho_default}")
    print(f"Number of passive elements read = {len(passive_ids)}")


if __name__ == "__main__":
    main()