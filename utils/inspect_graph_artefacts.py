import argparse
import os
import pandas as pd
from pathlib import Path

def display_parquet_summary(file_path: Path, file_description: str, head_n: int = 5):
    """Loads a Parquet file and displays its shape, head, and specific examples if applicable."""
    if not file_path.exists():
        print(f"\n--- {file_description} ({file_path.name}) ---")
        print(f"File not found: {file_path}")
        return

    try:
        df = pd.read_parquet(file_path)
        print(f"\n--- {file_description} ({file_path.name}) ---")
        print(f"Available Columns: {df.columns.tolist()}")

        if file_description == "Entities":
            name_cols = ['title', 'name', 'human_readable_id', 'label']
            desc_col = 'description'
            actual_name_col = None
            for col in name_cols:
                if col in df.columns:
                    actual_name_col = col
                    break
            
            if actual_name_col and desc_col in df.columns:
                print(f"\nExample Entities (Top {head_n}):")
                for i, row in df.head(head_n).iterrows():
                    title = row[actual_name_col]
                    description = row[desc_col]
                    print(f"  {i+1}. Title: {title}")
                    print(f"     Description: {description}")
            elif actual_name_col:
                print(f"\nExample Entity Names/IDs (from '{actual_name_col}' column, top {head_n}, description column not found):")
                unique_names = df[actual_name_col].dropna().astype(str)
                for i, name in enumerate(unique_names.head(head_n)):
                    print(f"  {i+1}. {name}")
            else:
                print("\nCould not find a typical 'name' or 'title' column for entities. Showing first few rows of all columns:")
                print(df.head(head_n).to_string())

        elif file_description == "Relationships":
            source_col, target_col, rel_desc_col = None, None, None
            src_candidates = ['source', 'src', 'source_id']
            tgt_candidates = ['target', 'tgt', 'target_id']
            desc_candidates = ['description', 'type', 'label', 'relationship_type', 'edge_type']
            
            for col in src_candidates:
                if col in df.columns: source_col = col; break
            for col in tgt_candidates:
                if col in df.columns: target_col = col; break
            for col in desc_candidates:
                if col in df.columns: rel_desc_col = col; break

            if source_col and target_col and rel_desc_col and 'weight' in df.columns:
                print(f"\nExample Relationships (Source -> Description -> Target, Weight, top {head_n}):")
                for index, row in df.head(head_n).iterrows():
                    src_val = row[source_col]
                    tgt_val = row[target_col]
                    desc_val = row[rel_desc_col]
                    weight_val = row['weight']
                    print(f"  {index+1}. Source: {src_val}")
                    print(f"     Target: {tgt_val}")
                    print(f"     Description: {desc_val}")
                    print(f"     Weight: {weight_val}")
            elif source_col and target_col and rel_desc_col:
                print(f"\nExample Relationships (Source -> Description -> Target, top {head_n}, weight column not found):")
                for index, row in df.head(head_n).iterrows():
                    src_val = row[source_col]
                    tgt_val = row[target_col]
                    desc_val = row[rel_desc_col]
                    print(f"  {index+1}. '{src_val}' --[{desc_val}]--> '{tgt_val}'")
            else:
                print("\nCould not find typical 'source', 'target', 'description/type', and/or 'weight' columns. Showing first few rows of all columns:")
                print(df.head(head_n).to_string())

        elif file_description == "Community Reports" and 'summary' in df.columns:
            print(f"\nTop {head_n if head_n <=2 else 2} Community Report Summaries (from 'summary' column):")
            if not df.empty:
                for i in range(min(len(df), head_n if head_n <=2 else 2)):
                    print(f"\nReport {i+1}:")
                    print(df['summary'].iloc[i])
            else:
                print("No community reports found in the file.")
        
        elif file_description == "Communities Definition":
            pass

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
        help="Number of example items to display (default: 5)."
    )
    parser.add_argument(
        "--show", 
        type=str, 
        nargs='*',
        choices=['entities', 'relationships', 'reports', 'all'],
        default=['all'],
        help="Specify which artefact(s) to show: entities, relationships, reports, or all (default: all)."
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

    if show_all or 'reports' in args.show:
        reports_file = artefacts_path / "community_reports.parquet" 
        display_parquet_summary(reports_file, "Community Reports", args.head if args.head <=2 else 2)

if __name__ == "__main__":
    main() 