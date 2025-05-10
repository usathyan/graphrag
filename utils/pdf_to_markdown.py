import argparse
import os
from pathlib import Path
from markitdown import MarkItDown

def convert_pdf_to_markdown(input_path_str: str, output_path_str: str = None):
    """
    Converts a PDF file or all PDF files in a directory to Markdown.

    Args:
        input_path_str: Path to the input PDF file or directory containing PDF files.
        output_path_str: Path to the output directory. Required if input_path_str is a directory.
                         If input_path_str is a file, output is saved in the same directory
                         with .md extension, and this argument is ignored unless explicitly provided
                         for a different output directory.
    """
    input_path = Path(input_path_str)
    md_converter = MarkItDown() # Initialize with default settings

    if not input_path.exists():
        print(f"Error: Input path does not exist: {input_path}")
        return

    if input_path.is_file():
        if input_path.suffix.lower() != '.pdf':
            print(f"Error: Input file is not a PDF: {input_path}")
            return

        if output_path_str:
            output_dir = Path(output_path_str)
            if not output_dir.is_dir():
                 print(f"Error: Specified output path '{output_dir}' for a single file input must be an existing directory or not provided (to save in same dir).")
                 return
            output_file_path = output_dir / f"{input_path.stem}.md"
        else:
            output_file_path = input_path.with_suffix('.md')
        
        output_file_path.parent.mkdir(parents=True, exist_ok=True)

        print(f"Converting file: {input_path} to {output_file_path}...")
        try:
            result = md_converter.convert(str(input_path))
            with open(output_file_path, 'w', encoding='utf-8') as f:
                f.write(result.text_content)
            print(f"Successfully converted: {output_file_path}")
        except Exception as e:
            print(f"Error converting file {input_path}: {e}")

    elif input_path.is_dir():
        if not output_path_str:
            print("Error: Output directory must be specified when input is a directory.")
            return

        output_dir = Path(output_path_str)
        output_dir.mkdir(parents=True, exist_ok=True)
        print(f"Processing PDF files in directory: {input_path}")
        print(f"Outputting Markdown to directory: {output_dir}")

        pdf_files_found = False
        for item in input_path.iterdir():
            if item.is_file() and item.suffix.lower() == '.pdf':
                pdf_files_found = True
                output_file_path = output_dir / f"{item.stem}.md"
                print(f"Converting file: {item} to {output_file_path}...")
                try:
                    result = md_converter.convert(str(item))
                    with open(output_file_path, 'w', encoding='utf-8') as f:
                        f.write(result.text_content)
                    print(f"Successfully converted: {output_file_path}")
                except Exception as e:
                    print(f"Error converting file {item}: {e}")
        
        if not pdf_files_found:
            print(f"No PDF files found in directory: {input_path}")
    else:
        print(f"Error: Input path is not a valid file or directory: {input_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert PDF files to Markdown using MarkItDown.")
    parser.add_argument("input_path", help="Path to the input PDF file or directory.")
    parser.add_argument("-o", "--output_path", help="Path to the output directory (required if input is a directory, optional for single file).", default=None)
    
    args = parser.parse_args()
    convert_pdf_to_markdown(args.input_path, args.output_path) 