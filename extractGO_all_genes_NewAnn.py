#!/usr/bin/env python3
"""extractGo_NewAnn.py

Stream a large tab-delimited file and extract Gene Ontology (GO) IDs.

Features:
- Read plain or gzipped files.
- Select column by index or name (supports header row).
- Regex-based GO ID extraction (default: GO:\\d+).
- Output modes: per-occurrence, unique, or counts.
- Memory-efficient streaming; counts kept in memory only when requested.

Example:
  python3 extractGo_NewAnn.py data.tsv --column go_terms --header --counts

"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from collections import Counter
from typing import Iterable, IO, Optional


GO_RE = re.compile(r"GO:\d+")


def open_maybe_gz(path: str) -> IO[str]:
	if path == "-":
		return sys.stdin
	if path.endswith(".gz"):
		return gzip.open(path, "rt", encoding="utf-8", errors="replace")
	return open(path, "r", encoding="utf-8", errors="replace")


def extract_from_field(field: str, pattern: re.Pattern) -> list[str]:
	return pattern.findall(field)


def iter_fields_from_lines(f: IO[str], delim: str, column: Optional[int], all_columns: bool, has_header: bool) -> Iterable[tuple[int, list[str]]]:
	"""Yield (line_number, fields) for each data line.

	If has_header is True, the first line is assumed header and skipped.
	"""
	lineno = 0
	first = True
	for raw in f:
		lineno += 1
		line = raw.rstrip("\n")
		if first and has_header:
			first = False
			header_fields = line.split(delim)
			continue
		first = False
		fields = line.split(delim)
		if not all_columns and column is not None:
			# safe get
			try:
				fields = [fields[column]]
			except Exception:
				fields = [""]
		yield lineno, fields


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
	p = argparse.ArgumentParser(description="Extract GO IDs from a large tab-delimited file (streaming)")
	p.add_argument("input", help="Input file path, or - for stdin")
	p.add_argument("--column", help="Column index (0-based) or name to search for GO IDs. If omitted, all columns are searched.")
	p.add_argument("--header", action="store_true", help="Treat first line as header (required if --column is a name).")
	p.add_argument("--delimiter", default="\t", help="Field delimiter (default: TAB)")
	p.add_argument("--pattern", default=r"GO:\d+", help="Regex pattern to match GO IDs (default: 'GO:\d+')")
	p.add_argument("--unique", action="store_true", help="Output unique GO IDs (one per line)")
	p.add_argument("--counts", action="store_true", help="Output GO ID counts as TSV: GO<TAB>count (implies --unique)")
	p.add_argument("--out", help="Output file (default stdout)")
	p.add_argument("--progress", type=int, default=100000, help="Print progress every N lines (default 100000). 0 to disable")
	p.add_argument("--all-columns", action="store_true", help="Search all columns instead of a single column")
	p.add_argument("--group-by-gene", action="store_true", help="Produce a two-column TSV: Gene<TAB>comma-separated-GO-IDs (unique)")
	p.add_argument("--gene-column", help="Column index (0-based) or name for gene identifier (default: 'gene' or first column)")
	p.add_argument("--go-columns", help="Comma-separated list of column names or 0-based indices to extract GO IDs from. If omitted, auto-detect columns whose header contains 'go' (case-insensitive).")
	return p.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> int:
	args = parse_args(argv)

	# determine column index if provided as integer
	col_index: Optional[int] = None
	col_name: Optional[str] = None
	if args.column is not None:
		try:
			col_index = int(args.column)
		except ValueError:
			col_name = args.column

	pattern = re.compile(args.pattern)

	out_f: IO[str]
	if args.out:
		out_f = open(args.out, "w", encoding="utf-8")
	else:
		out_f = sys.stdout

	counter: Optional[Counter[str]] = Counter() if args.counts else None

	with open_maybe_gz(args.input) as inf:
		header_fields: Optional[list[str]] = None
		# read header if requested
		if args.header:
			header_line = inf.readline()
			if not header_line:
				print("ERROR: Input appears empty while expecting header", file=sys.stderr)
				return 2
			header_fields = header_line.rstrip("\n").split(args.delimiter)

		# Resolve gene column
		gene_index: Optional[int] = None
		if args.gene_column is not None:
			try:
				gene_index = int(args.gene_column)
			except ValueError:
				if not header_fields:
					print("ERROR: --gene-column provided as name but --header not set", file=sys.stderr)
					return 4
				try:
					gene_index = header_fields.index(args.gene_column)
				except ValueError:
					print(f"ERROR: Gene column name '{args.gene_column}' not found in header", file=sys.stderr)
					return 3
		else:
			# default: if header exists try to find 'gene', else use column 0
			if header_fields:
				low = [h.lower() for h in header_fields]
				if 'gene' in low:
					gene_index = low.index('gene')
				else:
					gene_index = 0
			else:
				gene_index = 0

		# Resolve go columns
		go_col_indices: list[int] = []
		if args.go_columns:
			for token in args.go_columns.split(','):
				token = token.strip()
				if not token:
					continue
				try:
					go_col_indices.append(int(token))
				except ValueError:
					if not header_fields:
						print("ERROR: --go-columns provided as names but --header not set", file=sys.stderr)
						return 5
					try:
						go_col_indices.append(header_fields.index(token))
					except ValueError:
						print(f"ERROR: GO column name '{token}' not found in header", file=sys.stderr)
						return 6
		else:
			# auto-detect columns whose header contains 'go' (case-insensitive)
			if header_fields:
				for i, h in enumerate(header_fields):
					if 'go' in h.lower():
						go_col_indices.append(i)

		# If group mode and no GO columns found, warn and continue (will produce empty GO lists)
		if args.group_by_gene and not go_col_indices:
			print("WARNING: No GO columns auto-detected; output may be empty. Use --go-columns to specify columns.", file=sys.stderr)

		# iterate lines
		processed = 0
		# remaining lines: use iterator that does not skip header
		for lineno, fields in iter_fields_from_lines(inf, args.delimiter, col_index, args.all_columns, False):
			processed += 1
			if args.progress and processed % args.progress == 0:
				print(f"Processed {processed} lines...", file=sys.stderr)

			# Group-by-gene output
			if args.group_by_gene:
				# get gene id
				gene = ''
				try:
					gene = fields[gene_index]
				except Exception:
					gene = ''

				gos_set: set[str] = set()
				# if user specified go columns, use them; otherwise search all fields
				if go_col_indices:
					for gi in go_col_indices:
						try:
							fld = fields[gi]
						except Exception:
							fld = ''
						if not fld:
							continue
						for m in pattern.findall(fld):
							gos_set.add(m)
				else:
					# search all fields
					for fld in fields:
						if not fld:
							continue
						for m in pattern.findall(fld):
							gos_set.add(m)

				if gos_set:
					gos_list = sorted(gos_set)
					out_f.write(f"{gene}\t{','.join(gos_list)}\n")
				else:
					out_f.write(f"{gene}\t\n")
				continue

			# default per-field extraction behavior
			for fld in fields:
				if not fld:
					continue
				matches = pattern.findall(fld)
				if not matches:
					continue
				for m in matches:
					if counter is not None:
						counter[m] += 1
					else:
						# immediate write occurrence: one per line
						out_f.write(f"{m}\n")

	# if counts requested, write them sorted by count desc
	if counter is not None:
		# If user asked for unique only (counts implies unique output format)
		for go_id, cnt in counter.most_common():
			out_f.write(f"{go_id}\t{cnt}\n")

	if args.out:
		out_f.close()

	return 0


if __name__ == "__main__":
	raise SystemExit(main())
