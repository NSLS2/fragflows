import string
import pandas as pd
import pytz
import yaml
from pathlib import Path
import os

def letter_generator():
    def get_letter(num):
        result = []
        while num > 0:
            num, remainder = divmod(num - 1, 26)
            result.append(string.ascii_uppercase[remainder])
        return "".join(reversed(result))

    for k in range(1, 15000):
        yield get_letter(k)

def flatten_dict(d, parent_key="", sep="."):
    items = []
    for k, v in d.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k
        if isinstance(v, dict):
            items.extend(flatten_dict(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)

def convert_iso_date_to_ymd(date_str: str, target_tz: str) -> str:
    """Convert an ISO date string to a date string in the format YYYY-MM-DD, adjusting for the specified timezone."""

    dt = pd.to_datetime(date_str)  # Validate the input date string

    if target_tz:
        dt = dt.tz_convert(pytz.timezone(target_tz))

    return dt.strftime("%Y-%m-%d")


class PathResolver:
    def __init__(self, mapping_config: str = None):
        self.mappings = []
        self._load_mappings(mapping_config)

    def _load_mappings(self, mapping_config: str) -> dict:
        
        if mapping_config is None:
            mapping_config = "config.yaml"
        
        with open(mapping_config,"r") as f:
            config = yaml.safe_load(f)
            path_mappings = config.get("path_mappings", [])
            
            if len(path_mappings) == 0:
                return
            
            if len(path_mappings) % 2 != 0:
                raise ValueError("path_mappings should contain pairs of old and new root paths")
            
            for i in range(0, len(path_mappings), 2):
                old_root_path = path_mappings[i]
                new_root_path = path_mappings[i + 1]
                self.mappings.append({'old_root': old_root_path, 'new_root': new_root_path})

    def resolve(self, file_path: str) -> str:
        
        normalized_path = os.path.normpath(file_path)

        for mapping in self.mappings:
            old_root = os.path.normpath(mapping['old_root'])
            new_root = os.path.normpath(mapping['new_root'])
            
            if normalized_path.startswith(old_root):
                relative_path = os.path.relpath(normalized_path, old_root)
                resolved_path = os.path.join(new_root, relative_path)
                return resolved_path
        
        return file_path  # Return original if no mapping applies
            
    
    def resolve_dataframe_column(self, df, column: str, inplace: bool = True) -> None:
        """
        Apply path resolution to all entries in a pandas DataFrame column.
        
        Args:
            df: pandas DataFrame
            column: Column name containing filepaths
            inplace: If True, modify df in place; if False, return new df
        
        Returns:
            None if inplace=True, otherwise modified DataFrame
        """
        if inplace:
            df[column] = df[column].apply(self.resolve)
        else:
            return df.copy().assign(**{column: df[column].apply(self.resolve)})