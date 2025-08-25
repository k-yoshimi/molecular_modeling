from .build_system_from_toml import build_system_from_toml
from .write_structure import write_structure

def main():
    system = build_system_from_toml()
    write_structure(system)  
    return

if __name__ == "__main__":
    main()