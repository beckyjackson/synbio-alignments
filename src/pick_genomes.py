import csv
import os

from argparse import ArgumentParser


def main():
	parser = ArgumentParser()
	parser.add_argument("ani_directory")
	parser.add_argument("output")
	args = parser.parse_args()

	output = []
	for f in os.listdir(args.ani_directory):
		if not f.endswith(".tsv"):
			continue
		with open(os.path.join(args.ani_directory, f), "r") as fr:
			reader = csv.DictReader(fr, delimiter="\t")
			if "unicycler_assembly" not in reader.fieldnames:
				continue
			scores = {}
			for row in reader:
				genome = row[""]
				if "_assembly" in genome:
					# We only care about reference genomes
					continue
				scores[genome] = row["unicycler_assembly"]
			if not scores:
				continue
			scores = {k: v for k, v in sorted(scores.items(), key=lambda item: item[1], reverse=True)}
			output.append({"sample": f.split(".")[0], "genome": list(scores.keys())[0]})

	with open(args.output, "w") as fw:
		writer = csv.DictWriter(fw, lineterminator="\n", fieldnames=["sample", "genome"])
		writer.writeheader()
		writer.writerows(output)


if __name__ == '__main__':
	main()
