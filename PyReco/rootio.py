import uproot as root

def open_tree(filename: str, treename: str):
    file = root.open(filename)
    print(f"- \"{treename}\" from {filename} opened with success...")
    return file[treename]