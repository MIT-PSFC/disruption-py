from typing import Dict

def map_string_attributes_to_enum(obj, enum_translations : Dict):
    for field_name, enum_class in enum_translations.items():
            if hasattr(obj, field_name) and isinstance(getattr(obj, field_name), str):
                try:
                    enum_value = enum_class(getattr(obj, field_name))
                    setattr(obj, field_name, enum_value)
                except ValueError:
                    raise ValueError(f"Invalid enum value for field '{field_name}': {getattr(obj, field_name)}")