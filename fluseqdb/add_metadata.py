import json
import argparse
import pandas as pd
from .fluseqdb import FluSeqDatabase


def dicts_inconsistent(a: dict, b: dict) -> bool:
    shared_keys = set(a) & set(b)
    return any(a[k] != b[k] for k in shared_keys)


def add_metadata(db_dir: str, df: pd.DataFrame) -> None:
    fsdb = FluSeqDatabase(db_dir)

    for _, incoming_data in df.iterrows():

        for segment in fsdb.segments:

            isolate_id = incoming_data["isolate_id"]

            if not fsdb.exists(isolate_id=isolate_id, segment=segment):
                continue

            json_path = fsdb.path("sequences", segment, f"{isolate_id}.json")

            # If there is no existing metadata the output is the current row in the table
            if not json_path.exists():
                output = dict(incoming_data)

            # If there is existing metadata, check it is consistent, merge it and then output
            else:

                with open(json_path, "r") as fobj:
                    try:
                        existing_data = json.load(fobj)
                    except json.decoder.JSONDecodeError as err:
                        print(isolate_id, segment)
                        raise err

                if dicts_inconsistent(existing_data, incoming_data):
                    raise ValueError(
                        f"inconsistent data for {isolate_id}:\n"
                        f"existing data: {existing_data}\n"
                        f"incoming data: {incoming_data}"
                    )

                output = {**existing_data, **incoming_data}

            with open(json_path, "w") as fobj:
                json.dump(output, fobj)


def main():
    parser = argparse.ArgumentParser(
        "fsdb_addmeta", description="Add metadata to a fluseq db."
    )
    parser.add_argument("db_dir", help="Root directory of a fluseq db.")
    parser.add_argument(
        "table",
        help="CSV file containing metadata to add. First column contains isolate_id. Remaining "
        "columns are added to the isolate_id's data.",
    )
    args = parser.parse_args()

    add_metadata(db_dir=args.db_dir, df=pd.read_csv(args.table))
