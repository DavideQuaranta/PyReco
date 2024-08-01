import yaml

def read_configuration_file(name: str) -> dict[str, str]:
    """
    Reads a configuration file in YAML format and returns its contents as a dictionary.

    Args:
        name (str): The name or path of the configuration file to read.

    Returns:
        dict[str, str]: A dictionary containing the configuration key-value pairs.

    Raises:
        FileNotFoundError: If the specified file does not exist.
        yaml.YAMLError: If there is an error in parsing the YAML file.
        
    Example:
        >>> config = read_configuration_file('config.yaml')
        >>> print(config)
        {'key1': 'value1', 'key2': 'value2'}
    """
    with open(name, 'r') as file:
        config = yaml.safe_load(file)
    
    print(f"- configuration file \"{name}\" loaded...")

    return config



if __name__ == "__main__":
    config = read_configuration_file("../config.yaml")
    print(config)
