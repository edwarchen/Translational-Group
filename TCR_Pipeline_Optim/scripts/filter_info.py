#!/usr/bin/env python3
import argparse
import csv


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--sed-id", required=True)
    parser.add_argument("--adapter-pattern", required=True)
    parser.add_argument("--project-keyword", required=True)
    parser.add_argument("--chain", required=True)
    parser.add_argument("--receptor", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    with open(args.input, newline="", encoding="utf-8") as input_fh, open(
        args.output, "w", newline="", encoding="utf-8"
    ) as output_fh:
        reader = csv.reader(input_fh)
        writer = csv.writer(output_fh)
        writer.writerow(
            [
                "Sed_ID",
                "Sample_ID",
                "Lib_number",
                "Name",
                "Panel",
                "Data_ID",
                "Chain",
                "Receptor",
            ]
        )
        for row in reader:
            if (
                len(row) >= 24
                and row[0] == args.sed_id
                and args.adapter_pattern in row[11]
                and args.project_keyword in row[7]
            ):
                writer.writerow(
                    [row[0], row[1], row[2], row[3], row[15], row[23], args.chain, args.receptor]
                )


if __name__ == "__main__":
    main()
