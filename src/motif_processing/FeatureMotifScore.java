package motif_processing;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.HashMap;

import jbfunctions.GenomeFunctions;
import jbfunctions.PWMFunctions;

public class FeatureMotifScore {

	public static void main(String[] args) throws IOException {
		
		// Set the input file name
		String input_file_name = args[0];
		
		// Set the output file name
		String output_file_name = args[1];
		
		// Set the input ACS pwm file
		String input_pwm_file_name = args[2];
		
		// Set the bg file name
		String bg_model_file_name = args[3];
				
		// Set the genome fasta header
		String genome_fasta_header = args[5];
		
		////////////////////////////////////////////////////////////////////////////////////////////////
				
		// Open the output buffer
		BufferedWriter output = new BufferedWriter(new FileWriter(output_file_name));
		
		// Write the header
		output.write("name,chr,left_coord,right_coord,strand,log2_motif_score,seq\n");
		
		// Read in the acs pwm
		double[][] tf_pwm = PWMFunctions.read_in_pwm(input_pwm_file_name);
		
		// Get the background model
		HashMap<String, Double> bg_model = PWMFunctions.read_in_bg_model(bg_model_file_name);
			
		// Set the DecimalFormat
		DecimalFormat df = new DecimalFormat("#.####");
		
		// Open the input buffer
		BufferedReader input = new BufferedReader(new FileReader(input_file_name));
		
		// Initialize the chromosome
		String chr = "1";
		
		String chr_seq =  GenomeFunctions.read_in_chr_fasta_file(chr, genome_fasta_header);
		
		// Read in the header
		String line = input.readLine();
		
		while((line = input.readLine()) != null){
			
			// Split the feature
			String[] line_arr = line.split(",");
			
			// Get the chr
			String feature_chr = line_arr[1];
						
			// Read in the chr_seq if necessary
			if(!feature_chr.equals(chr)){
				chr = feature_chr;
				chr_seq =  GenomeFunctions.read_in_chr_fasta_file(feature_chr, genome_fasta_header);
			}
						
			// Get the genome coordinates
			int left_coord = Integer.parseInt(line_arr[2]);
			int right_coord = Integer.parseInt(line_arr[3]);
			
			// Get the strand
			char str = line_arr[4].charAt(0);
			
			int motif_pos = left_coord;
			if(str == '-'){
				motif_pos = right_coord;
			}
			
			// Find the log_diff_score
			double log_diff_score = PWMFunctions.get_log_score(chr_seq, tf_pwm, motif_pos, str, bg_model);
		
			// Write the output
			output.write(line + "," + 
						 df.format(log_diff_score) + "," +
						 PWMFunctions.get_seq_to_evaluate(chr_seq, left_coord, right_coord, str) + "\n"
						);
				
								
		}
		
		// Close the input buffer
		input.close();
		
		// Close the output buffer
		output.close();
		
	}

}
