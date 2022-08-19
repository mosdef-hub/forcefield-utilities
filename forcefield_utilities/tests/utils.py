from pathlib import Path


def get_test_file_path(filename):
    """Give a filename, return its location in test files."""
    file_path = Path(__file__).parent / "files" / filename

    if not file_path.resolve().exists():
        raise FileNotFoundError(
            f"File {filename} not found in {file_path.parent}"
        )

    return str(file_path)
