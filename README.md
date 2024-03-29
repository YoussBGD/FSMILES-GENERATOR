Introduction:
-

This Python script reads 3D molecular structures either from local SDF files or directly via PubChem CID, and gives its FSMILES representation.



Dependencies:
-
The required packages are:

  - RDKit: A collection of cheminformatics and machine learning tools.

  - NetworkX: A package for the creation, manipulation, and study of the structure, dynamics, and functions of complex networks.

  - Requests: A simple HTTP library for Python.

  - Requests-cache: A transparent persistent cache for the Requests library.

You can install these dependencies using pip  by running the following commands in your terminal:
  pip install rdkit networkx requests requests-cache
 


Input Data:
-
The script can process molecules from two sources:

  - Local SDF files containing the molecular structure.
  - PubChem Compound ID (CID).

Run the Script: 
-
Navigate to the directory containing the script and run it using Python:

  python FSMILES_representation.py "input sdf file or pubchem CID" "output_file.txt"



Output:
-
The script will generate a .txt file  containing the FSMILES representation of your molecule.


What the Code Does:
-

- Reads Molecules: From SDF files or via PubChem CIDs, converting them to RDKit molecule objects.

- Identifies Cuttable Bonds: Finds bonds in the molecule that can be cut based on specific criteria (not in rings, not hydrogen atoms, etc.).

- Fragments Molecules: Fragments the molecule at the identified cuttable bonds, replacing them with dummy atoms to indicate points of fragmentation.
- Constructs a Connectivity Graph: Builds a graph representing the connectivity of the fragments.
- Performs DFS Traversal: Traverses the graph using Depth-First Search (DFS) starting from a root atom to maintain the structural integrity of the molecule in the SMILES representation.
- Reconstructs the modified SMILE: Creates a modified SMILES string (FSMILES) based on the traversal, annotating cycle sizes for each atom.
- Outputs Modified SMILES String: Saves the  FSMILES string to your output.txt file.


