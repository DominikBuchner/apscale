# main function of the swarm clustering script
def main(project=Path.cwd()):
    """Main function of the script. Default values can be changed via the input file.
    Will perform swarm clustering on the individual files. If run without denoising will also perform chimera removal and
    exchange the fasta headers with sha 256 hashes for easier processing.

    Args:
        project (str, optional): Path to the apscale project. Defaults to Path.cwd().
    """
    print("Hello")
