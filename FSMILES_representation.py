from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem import rdmolfiles
import networkx as nx
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG
import re
import requests
import requests_cache
import sys
# Enable caching for requests to reduce the number of network calls for repeated molecules
requests_cache.install_cache('pubchem_cache', expire_after=86400)

#-----------------------------------------All_functions----------------------------------------


def read_mol(sdf_path):

    if any(char.isalpha() for char in sdf_path) or '.sdf' in sdf_path:
        try:
            supplier = Chem.SDMolSupplier(sdf_path)
            mol3D = next(iter(supplier), None)  # Get the first molecule from the supplier if more than one molecule.
            if mol3D is None:
                raise ValueError("\n\nNo molecules found in the file.\n\n")
                exit()
            smiles = Chem.MolToSmiles(mol3D)
        except Exception as e:
            print("\n\nInvalid SDF path or corrupted file.\n\n")
            exit()

        
        smiles = Chem.MolToSmiles(mol3D)  # Convert 3D mol to SMILES string.
    else :
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{sdf_path}/SDF"
    
        # Send a request to PubChem to get the molecule's SDF.

        try:
            response = requests.get(url)
            response.raise_for_status()  # This will raise an HTTPError for bad responses
            # Proceed with handling the response if no error is raised
        except requests.exceptions.HTTPError as e:
            print(f"\n\nFailed to retrieve data from PubChem: {e}\n\n")
            exit()
        # Use RDKit to read the SDF from the response text.
        sdf = Chem.MolFromMolBlock(response.text)

        # Convert the molecule to a SMILES string.
        smiles = Chem.MolToSmiles(sdf)
        
    mol = Chem.MolFromSmiles(smiles)  # Convert SMILES string back to mol for processing.
    return mol

    

def find_cuttable_bonds(mol):
    # Initialize lists to hold indices of cuttable bonds and the bond objects themselves.
    cuttable_bondsIDX = []  # To store indices of all cuttable bonds.
    cuttable_bonds = []  # To store the actual bond objects that can be cut.

    # Iterate over all bonds in the molecule to determine if they can be cut.
    for bond in mol.GetBonds():

        if (bond.GetBondType() == Chem.rdchem.BondType.SINGLE and  # Bond is a single bond.
            not bond.IsInRing() and  # Bond is not part of a ring.
            bond.GetBeginAtom().GetAtomicNum() != 1 and  # The starting atom is not hydrogen.
            bond.GetEndAtom().GetAtomicNum() != 1 and  # The ending atom is not hydrogen.
            (mol.GetRingInfo().NumAtomRings(bond.GetBeginAtomIdx()) > 0 or  # The starting atom is part of a ring.
             mol.GetRingInfo().NumAtomRings(bond.GetEndAtomIdx()) > 0)):  # The ending atom is part of a ring.
            cuttable_bondsIDX.append(bond.GetIdx())  # Append bond index to list of cuttable bond indices.
            cuttable_bonds.append(bond)  # Append the bond object to the list of cuttable bonds.

    # After identifying all cuttable bonds, extract the atom indices for each bond.
    cutbonds_atom_idx = []
    for bond in cuttable_bonds:
        # For each cuttable bond, store the indices of the begin and end atoms.
        cutbonds_atom_idx.append([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])
    
    # Return both the indices of the cuttable bonds and the atom indices for each cuttable bond.
    return cuttable_bondsIDX, cutbonds_atom_idx


def fragment_molecule_and_get_info(mol, bond_indices):
    """
    Fragment the input molecule at specified bond indices and extract information for each fragment.
    
    Parameters:
    - mol: The RDKit molecule to be fragmented.
    - bond_indices: A list of bond indices at which the molecule will be fragmented.
    
    Returns:
    - frag_atoms: A dictionary where each key is a fragment identifier (e.g., 'frag1') and the value
      is a tuple containing the atom indices list and the SMILES string of the fragment.
    """
    # Perform fragmentation
    if bond_indices:
        fragments = rdmolops.FragmentOnBonds(mol, bond_indices, addDummies=True)
    else:
        fragments = mol #if there in no cuttable bonds (no cycles) ex : cid 545889

    frag_mols = Chem.GetMolFrags(fragments, asMols=True, sanitizeFrags=True) # A list of RDKit molecule objects for each fragment.
    
    # Extract fragment information
    atom_indices_in_fragments = Chem.GetMolFrags(fragments, asMols=False, sanitizeFrags=False)
    frag_atoms = {}
    for idx, (frag_mol, atom_indices) in enumerate(zip(frag_mols, atom_indices_in_fragments), start=1):
        frag_atoms[f'frag{idx}'] = (atom_indices, Chem.MolToSmiles(frag_mol, isomericSmiles=True))
    
    return frag_atoms


def get_atom_cycle_info(mol):
    """
    Determine the cycle size for each atom in the molecule.
    
    Iterates over all atoms in the molecule and checks if each atom is part of a cycle.
    If so, it determines the sizes of cycles that the atom belongs to.
    The function then stores and returns this information in a dictionary where:
    - The key is a string combining the atom's symbol and its index.
    - The value is another dictionary with two keys: 'cycle_size' (size of the first cycle the atom belongs to or 0 if not in a cycle) and 'index' (atom's index in the molecule).
    """
    atom_info = {}
    for atom in mol.GetAtoms():
        atom_idx = atom.GetIdx()
        in_cycle = atom.IsInRing()
        cycle_sizes = [len(cycle) for cycle in mol.GetRingInfo().AtomRings() if atom_idx in cycle]
        cycle_size = cycle_sizes[0] if cycle_sizes else 0
        atom_info[atom.GetSymbol() + str(atom_idx)] = {"cycle_size": cycle_size, "index": atom_idx}
    return atom_info



def find_fragment(atom_idx, frag_atoms):
    """
    Finds the fragment containing a given atom index.
    Args:
        atom_idx (int): The index of the atom.
        frag_atoms (dict): A dictionary containing information about fragments.
    Returns:
        frag_name (str): The name of the fragment containing the atom.
    """
    for frag_name, (atom_indices, _) in frag_atoms.items():
        if atom_idx in atom_indices:
            return frag_name
    return None

def find_root_atom(mol):
    """
    Identifies the root atom of the molecule, defined as the non-ring atom with the highest degree of connectivity, the DFS algo will start from the fragment containing this atom.
    Args:
        mol (Chem.Mol): The RDKit molecule object.
    Returns:
        root_atom_idx (int): The index of the root atom.
    """
    max_degree = 0
    root_atom_idx = None
    for atom in mol.GetAtoms():
        if mol.GetRingInfo().NumAtomRings(atom.GetIdx()) == 0:  # Check if atom is not in a ring.
            degree = len(atom.GetNeighbors())
            if degree > max_degree:
                max_degree = degree
                root_atom_idx = atom.GetIdx()
    return root_atom_idx

def find_fragment_key_with_root_atom(frag_atoms, root_atom_idx):
    """
    Finds the fragment key that contains the root atom.
    Args:
        frag_atoms (dict): A dictionary containing information about fragments.
        root_atom_idx (int): The index of the root atom.
    Returns:
        frag_key (str): The key of the fragment containing the root atom, if found.
    """
    for frag_key, (atom_indices, _) in frag_atoms.items():
        if root_atom_idx in atom_indices:
            return frag_key
    return None

def dfs(graph, start, visited=None, path=None):
    """
    Performs a Depth-First Search (DFS) on the graph starting from a given node (fragment containing root atom).
    Args:
        graph (nx.Graph): The graph to traverse.
        start (str): The starting node of the DFS.
        visited (set, optional): A set of visited nodes.
        path (list, optional): The path traversed so far.
    Returns:
        path (list): The complete path after DFS traversal.
    """
    if visited is None:
        visited = set()
    if path is None:
        path = []
    visited.add(start)
    path.append(start)
    for neighbor in graph.neighbors(start):
        if neighbor not in visited:
            dfs(graph, neighbor, visited, path)
    return path

def reconstruct_modified_dfs_smile(graph, start_fragment, frag_atoms):
    """
    Reconstructs the modified SMILES string from the graph using DFS traversal.
    Args:
        graph (nx.Graph): The graph representing fragment connectivity.
        start_fragment (str): The starting fragment for the DFS.
        frag_atoms (dict): A dictionary containing information about fragments.
    Returns:
        fsmile_cl (str): The reconstructed SMILES string with modifications.
    """
    visited_fragments = dfs(graph, start_fragment)
    fsmile = "'start'"
    for frag_key in visited_fragments:
        fragment_smile = frag_atoms[frag_key][1]  # Directly access the SMILE of the fragment.
        fsmile += fragment_smile + "'sep'"
    fsmile += "'end'"  # Append 'end' after traversing all fragments.
    fsmile_cl = re.sub(r"sep'\[\d+\*\]|sep'\*", "sep'", fsmile)  # Clean up the SMILE string.
    return fsmile_cl

def construct_graph(frag_atoms, cutbonds_idx):
    """
    Constructs a graph based on the fragments and their connectivity.
    Args:
        frag_atoms (dict): A dictionary containing information about fragments.
        cutbonds_idx (list): A list of atom index pairs representing cuttable bonds.
    Returns:
        G (nx.Graph): The constructed graph representing fragment connectivity.
    """
    G = nx.Graph()
    # Add nodes for each fragment.
    for frag_name in frag_atoms.keys():
        G.add_node(frag_name)
    # Add edges based on cuttable bonds to represent fragment connectivity.
    for start_idx, end_idx in cutbonds_idx:
        frag1 = find_fragment(start_idx, frag_atoms)
        frag2 = find_fragment(end_idx, frag_atoms)
        if frag1 and frag2:
            G.add_edge(frag1, frag2)
    return G

def move_and_transform_consecutive_dummy_atoms_correctly(smiles):
    # Use a regular expression to identify the first occurrence of consecutive dummy atoms ([x*], where x is a digit)
    # immediately following the 'start' marker in the SMILES string.
    pattern = r"(?<='start')(\[\d+\*\])+"
    match = re.search(pattern, smiles)
    
    # If such a pattern is found in the SMILES string (only for the part succeding start 
    #(because of some non atomic elements which can be found after start and damage the smile representation )): 
    if match:
        # Extract the entire sequence of dummy atoms ([x*]) and record its start and end positions.
        dummy_atoms_sequence = match.group()
        start_pos = match.start()
        end_pos = match.end()

        # Remove the original sequence of dummy atoms from the SMILES string.
        smiles = smiles[:start_pos] + smiles[end_pos:]

        # Transform each dummy atom in the extracted sequence from [x*] to ([*]),
        # effectively removing the numeric identifier to standardize the representation.
        transformed_sequence = re.sub(r'\[\d+\*\]', r'([*])', dummy_atoms_sequence)

        # Find the position of the first actual atom (denoted by an alphabetical character)
        # that comes after the 'start' marker in the modified SMILES string.
        first_atom_pos_after_start = re.search(r'[A-Za-z]', smiles[len("'start'"):]).start() + len("'start'")

        # Insert the transformed sequence of dummy atoms right after the first actual atom found.
        smiles = smiles[:first_atom_pos_after_start + 1] + transformed_sequence + smiles[first_atom_pos_after_start + 1:]
    
    # Clean up the SMILES string by replacing any non-atomic element immediately following 'start'
    # and replacing all remaining instances of [x*] with [*] to standardize dummy atom representation.
    smiles = re.sub(r"(start')[^A-Za-z]+", r"\1", smiles)  # Clean up non-atomic elements after 'start'.
    cleaned_smiles = re.sub(r'\[\d+\*\]', '[*]', smiles)  # Replace all [x*] with [*].

    # Return the processed SMILES string with correctly transformed and placed dummy atoms.
    return cleaned_smiles


def annotate_with_cycle_size(cleaned_elements, fragment_info):
    # Initialize a list to hold elements after annotation.
    annotated_elements = []
    # Initialize an atom index counter to keep track of atoms within each fragment.
    atom_index = -1
    # Start with the assumption that we're working on the first fragment.
    current_fragment = "Fragment 1"

    # Iterate over each element in the list of cleaned elements extracted from the modified SMILES string.
    for element in cleaned_elements:
        # Handle special markers ('start', 'sep', 'end') which don't increment the atom index.
        if element in ["'start'", "'sep'", "'end'"]:
            # Add a '_0' suffix to these special elements for correct annotation.
            special_annotation = element.rstrip("'") + "_0'"
            annotated_elements.append(special_annotation)
            # When encountering a 'sep' marker, it indicates moving to the next fragment.
            if element == "'sep'":
                # Calculate the next fragment number and reset the atom index for the new fragment.
                current_fragment_number = int(current_fragment.split()[-1]) + 1
                current_fragment = f"Fragment {current_fragment_number}"
                atom_index = -1
        # Handle dummy atom symbols ([*]) and ([*]) which do increment the atom index but don't have cycle sizes.
        elif element in ["([*])", "[*]"]:
            atom_index += 1
            # Annotate these with "_0" to indicate no cycle size.
            annotated_elements.append(element + "_0")
        # For atomic symbols, increment the atom index and retrieve the cycle size from the fragment information.
        elif element.isalpha() and element not in ["@","H","h"]:
            atom_index += 1
            # Generate a unique key for the atom based on its symbol and index.
            if len(element) < 2:
                atom_key = f"{element.upper()}{atom_index}"
            else:
                atom_key = f"{element}{atom_index}"
            # Retrieve the cycle size information for the atom, defaulting to 0 if not found.
            cycle_size_info = fragment_info.get(current_fragment, {}).get(atom_key, {"cycle_size": 0})
            cycle_size = cycle_size_info.get("cycle_size", 0)
            # Annotate the element with its cycle size.
            cycle_annotation = f"_{cycle_size}" if cycle_size > 0 else "_0"
            annotated_elements.append(element + cycle_annotation)
        else:
            # For any other elements (not atomic or special markers), they don't modify the index.
            annotated_elements.append(element + "_0")

    # Combine the annotated elements into a single string and return it.
    return ''.join(annotated_elements)



def atoms_in_mol(mol):
    all_atoms_in_mol = []
    # Iterate through each atom in the molecule
    for atom in mol.GetAtoms():
        all_atoms_in_mol.append(atom.GetSymbol())
    unique_atoms = set(all_atoms_in_mol)
    unique_atoms = [atom for atom in unique_atoms if len(atom) > 1] #keep only atoms having more than one letter 
    return unique_atoms


def process_molecule_from_sdf(sdf_path):
    # Read the molecule from the specified SDF file.
    mol = read_mol(sdf_path)
    atoms_more_than_one_caracter = atoms_in_mol(mol)
    # Find indices of bonds that can be cut, based on certain criteria (e.g., not in rings, not connecting hydrogen atoms, etc.).
    cuttable_bonds_idx, cuttable_bonds_atom_idx = find_cuttable_bonds(mol)
    # Fragment the molecule at these cuttable bonds, replacing them with dummy atoms to indicate the points of fragmentation.
    #Collect information about these fragments, including atom indices within each fragment and their SMILES representation.
    frag_atoms = fragment_molecule_and_get_info(mol, cuttable_bonds_idx)
    # Construct a graph that represents the connectivity of these fragments, using the information about cuttable bonds.
    G = construct_graph(frag_atoms, cuttable_bonds_atom_idx)
    # Identify the 'root' atom of the molecule, which is used as a starting point for certain operations like DFS traversal.
    root_atom_idx = find_root_atom(mol)
    # Return the graph representing fragment connectivity, information about the fragments, and the index of the root atom.
    return G, frag_atoms, root_atom_idx, atoms_more_than_one_caracter


def process_new_smile_reconstruction(G, frag_atoms, root_atom_idx):
    # Find the fragment that contains the root atom. This fragment serves as the starting point for reconstructing the SMILES string.
    frag_key_with_root = find_fragment_key_with_root_atom(frag_atoms, root_atom_idx)
    # Reconstruct the SMILES string by traversing the graph starting from the root fragment, thereby ensuring that the structure is represented in a linear, coherent manner.
    FSMILE1 = reconstruct_modified_dfs_smile(G, frag_key_with_root, frag_atoms)
    # Perform additional processing on the reconstructed SMILES string to correctly position and format dummy atoms.
    modified_smiles = move_and_transform_consecutive_dummy_atoms_correctly(FSMILE1)
    # Return the processed SMILES string.
    return modified_smiles


def clean_smiles_string(modified_smiles):
    # Remove the 'start' and 'end' markers from the SMILES string, which are used for indicating the beginning and end of the molecular structure.
    fsmiles_clean = modified_smiles.replace("'start'", "").replace("'end'", "")
    # Split the cleaned SMILES string into fragments based on the 'sep' marker, which separates different parts of the molecule, potentially resulting from fragmentation.
    fragments = fsmiles_clean.split("'sep'")
    # Ensure the list of fragments does not include an empty string at the end, which could happen if 'sep' is at the end of the SMILES.
    if fragments[-1] == "":
        fragments.pop()
    # Return the list of fragment SMILES strings.
    return fragments


def extract_mod_fragment_informations(fragments):
    # Initialize a dictionary to hold information about each fragment, keyed by a fragment identifier and containing cycle size information.
    fragments_info = {}
    for i, smiles in enumerate(fragments, start=1):
        # Convert each fragment's SMILES string back into a molecule object.
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # If conversion is successful, get information about atom cycles within the fragment.
            atom_info = get_atom_cycle_info(mol)
            # Store this information, indexed by the fragment number.
            fragments_info[f"Fragment {i}"] = atom_info
        else:
            # If conversion fails, report an error for this fragment.
            print(f"Error converting SMILES to molecule for fragment {i}")
    return fragments_info



def compile_cleaned_elements(modified_smiles, atoms_more_than_one_caracter):

    atoms_pattern_part = "|".join(atoms_more_than_one_caracter)
    #print(all_atoms)
    # Original pattern with the dynamic part included for matching atom types
    pattern = rf"'start'|'end'|'sep'|\(\[\*\]\)|\[\*\]|{atoms_pattern_part}|[^']*?" # Use a regular expression to identify special markers ('start', 'end', 'sep', '([*])', '[*]') and other characters within the SMILES string.

    # Extract all elements that match this pattern from the modified SMILES string.
    elements = re.findall(pattern, modified_smiles)
    #print(elements)
    # Filter out any elements that are empty strings, to clean up the list.
    return [element for element in elements if element.strip()]


def prepare_smile_to_fsmile(sdf_path):
    # Integrate various steps to prepare a modified SMILES string from an SDF file, including processing the molecule, reconstructing the SMILES string, and cleaning it.
    G, frag_atoms, root_atom_idx,atoms_more_than_one_caracter = process_molecule_from_sdf(sdf_path)
    modified_smiles = process_new_smile_reconstruction(G, frag_atoms, root_atom_idx)
    fragments = clean_smiles_string(modified_smiles)
    fragments_info = extract_mod_fragment_informations(fragments)
    cleaned_elements = compile_cleaned_elements(modified_smiles,atoms_more_than_one_caracter)
    #Return the cleaned elements and information about the fragments for the last processing of adding cycles size.
    return cleaned_elements, fragments_info



def main(sdf_path, output_file):
    cleaned_elements, fragments_info = prepare_smile_to_fsmile(sdf_path)
    Fsmile = annotate_with_cycle_size(cleaned_elements, fragments_info)
   
    with open(output_file, "w") as output:
        output.write("Fsmile representation for molecule : {} \n\n".format(sdf_path))
        output.write(Fsmile)
    print("\n\n\nFinish processing, find the FSMILES representation of your molecule in the file {}".format(output_file))

#-----------------------------------End_fuctions--------------------------------------------    

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python code.py <input_path> <output_file>")
        sys.exit(1)
   
    input_path = sys.argv[1]
    output_file = sys.argv[2]
    main(input_path.strip(), output_file.strip())
