import re

def parse_treefile(treefile_content):
    # Regular expression to extract the identifier
    pattern = re.compile(r'(\d+_0_\w+)')
    identifiers = pattern.findall(treefile_content)
    return identifiers

def write_fasta(identifiers, output_fasta):
    with open(output_fasta, 'w') as fasta_file:
        for identifier in identifiers:
            fasta_file.write(f">{identifier}___\n")

def transform_identifiers(identifiers):
    # Placeholder for any transformation you want to apply to identifiers
    transformed_identifiers = [id for id in identifiers]  # No transformation in this example
    return transformed_identifiers

# Main process
treefile = "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/ITPKC/ITPKC.treefile"  # Replace with the actual path to your treefile
output_fasta = "/home/pythagoras/Documents/PhD/Evolution/Proteins/ITPKA/ITPKC/ITPKC_ownAnnotation.fasta"

# Read treefile content
with open(treefile, 'r') as file:
    treefile_content = file.read()

# Parse the treefile
identifiers = parse_treefile(treefile_content)

# Transform identifiers if needed
transformed_identifiers = transform_identifiers(identifiers)

# Write to FASTA
write_fasta(transformed_identifiers, output_fasta)

print(f"FASTA file '{output_fasta}' generated successfully!")
