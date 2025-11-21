import argparse
from Bio import SeqIO
import csv
from pathlib import Path
from Bio.SeqRecord import SeqRecord

def custom_abi_trim(seq_record, cutoff=0.05):
    """
    Trims the sequence using Richard Mott's modified trimming algorithm.
    Adapted from Bio.SeqIO.AbiIO._abi_trim
    """
    start = False  # flag for starting position of trimmed sequence
    segment = 20   # minimum sequence length
    trim_start = 0 # init start index

    if len(seq_record) <= segment:
        return seq_record
    else:
        # calculate base score
        score_list = [
            cutoff - (10 ** (qual / -10.0))
            for qual in seq_record.letter_annotations["phred_quality"]
        ]

        # calculate cumulative score
        cummul_score = [0]
        for i in range(1, len(score_list)):
            score = cummul_score[-1] + score_list[i]
            if score < 0:
                cummul_score.append(0)
            else:
                cummul_score.append(score)
                if not start:
                    # trim_start = value when cumulative score is first > 0
                    trim_start = i
                    start = True

        # trim_finish = index of highest cumulative score
        trim_finish = cummul_score.index(max(cummul_score))

        return seq_record[trim_start:trim_finish]
    



def generate_trim_report(pre_trim:SeqRecord,
                        post_trim:SeqRecord,
                        report_path:str,
                        sample_name:str
                        ) -> None:
    """
    Calculates quality metrics for pre and post trim and writes them to a CSV file.
    """

    # Protects from ZeroDivisionError
    if len(post_trim) > 0:
        post_avg = round(sum(post_trim.letter_annotations["phred_quality"]) / len(post_trim), 3)
    else:
        post_avg = 0.0

    report_dict = {
        "sample": sample_name,

        "removed_bp": len(pre_trim) - len(post_trim),
        
        "pre_trim_avg_phred_score": round(
            sum(pre_trim.letter_annotations["phred_quality"]) / len(pre_trim), 3
        ),
        
        "pre_trim_q20_count": sum(
            q > 20 for q in pre_trim.letter_annotations["phred_quality"]
        ),

        "post_trim_avg_phred_score": post_avg,
        
        "post_trim_q20_count": sum(
            q > 20 for q in post_trim.letter_annotations["phred_quality"]
        )
    }

    with open(report_path, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=report_dict.keys())
        writer.writeheader()
        writer.writerow(report_dict)
    


def main():
    # Parse arguments for Command-Line integration
    parser = argparse.ArgumentParser(description="Trim one .ab1 file using Mott's algorithm, \
                                     generating a .fastq trimmed read and a .csv report.")
    parser.add_argument("-i", "--input", 
                        required=True,
                        help="The path to the .ab1 file to be trimmed")
    parser.add_argument("-o", "--output_fastq", 
                        required=True,
                        help="What the generated trimmed file will be named. Must include the extension ('.fastq')")
    parser.add_argument("--output_report", 
                        help="The filename for the report file. Must include the extension ('.csv')")
    parser.add_argument("--trim_cutoff",
                        help="The cutoff value used by Mott's algorithm. Standard is 0.05, which corresponds to Q13.")
    args = parser.parse_args()

    # Read the ab1 file
    record_pre_trim = SeqIO.read(args.input, "abi")
    # Generate the trimmed record
    record_post_trim = custom_abi_trim(record_pre_trim, cutoff=float(args.trim_cutoff))

    SeqIO.write(record_post_trim, args.output_fastq, "fastq")

    input_path = Path(args.input)

    generate_trim_report(
        pre_trim=record_pre_trim,
        post_trim=record_post_trim,
        report_path=args.output_report,
        sample_name=input_path.stem
    )


if __name__ == "__main__":
    main()