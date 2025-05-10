import argparse
import os
import pandas as pd
from pathlib import Path

def display_parquet_summary(file_path: Path, file_description: str, head_n: int = 5):
    """Loads a Parquet file and displays its shape and head."""
    if not file_path.exists():
        print(f"\n--- {file_description} ({file_path.name}) ---")
        print(f"File not found: {file_path}")
        return

    try:
        df = pd.read_parquet(file_path)
        print(f"\n--- {file_description} ({file_path.name}) ---")
        print(f"Shape: {df.shape}")
        print(f"Head:\n{df.head(head_n)}")
        if file_description == "Community Reports" and 'summary' in df.columns:
            print("\nFirst Community Report Summary:")
            if not df.empty:
                print(df['summary'].iloc[0])
            else:
                print("No community reports found in the file.")

    except Exception as e:
        print(f"\n--- {file_description} ({file_path.name}) ---")
        print(f"Error reading or displaying file {file_path}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Inspect GraphRAG Parquet artefacts.")
    parser.add_argument(
        "artefacts_dir", 
        type=str, 
        help="Path to the GraphRAG output artefacts directory (e.g., output/YYYYMMDD-HHMMSS/artefacts/)."
    )
    parser.add_argument(
        "--head", 
        type=int, 
        default=5, 
        help="Number of rows to display from the head of each file (default: 5)."
    )
    parser.add_argument(
        "--show", 
        type=str, 
        nargs='*',
        choices=['entities', 'relationships', 'communities', 'reports', 'all'],
        default=['all'],
        help="Specify which artefact(s) to show: entities, relationships, communities, reports, or all (default: all)."
    )

    args = parser.parse_args()

    artefacts_path = Path(args.artefacts_dir)

    if not artefacts_path.is_dir():
        print(f"Error: Artefacts directory not found: {artefacts_path}")
        return

    print(f"Inspecting artefacts in: {artefacts_path}\n")

    show_all = 'all' in args.show

    if show_all or 'entities' in args.show:
        entities_file = artefacts_path / "entities.parquet"
        display_parquet_summary(entities_file, "Entities", args.head)

    if show_all or 'relationships' in args.show:
        relationships_file = artefacts_path / "relationships.parquet"
        display_parquet_summary(relationships_file, "Relationships", args.head)

    if show_all or 'communities' in args.show:
        communities_file = artefacts_path / "communities.parquet"
        display_parquet_summary(communities_file, "Communities Definition", args.head)
    
    if show_all or 'reports' in args.show:
        reports_file = artefacts_path / "community_reports.parquet"
        display_parquet_summary(reports_file, "Community Reports", args.head)

if __name__ == "__main__":
    main() 