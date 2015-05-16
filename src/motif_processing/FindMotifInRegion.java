package Motifs;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import jbfunctions.GenomeFunctions;
import jbfunctions.PWMFunctions;



public class FindMotifInRegion {

	public static ArrayList<FoundMotif> get_motif_list(
		String name, String chr, int start, int end, 
		double[][] pwm, String chr_seq, HashMap<String, Double> bg_model, double thresh_score
		){
		
		// Set the storage ArrayList
		ArrayList<FoundMotif> motif_list = new ArrayList<FoundMotif>();
		
		// Get the str array
		char[] str_arr = {'+', '-'};
					
		// Iterate over each strand
		for(char str : str_arr){
							
			for(int p = Math.max(start, pwm[0].length); 
					p <= Math.min(end, chr_seq.length() - pwm[0].length);
					p++
				){ 
					
				// Calculate the log-odds score
				double log_diff_score = PWMFunctions.get_log_score(
					chr_seq, pwm, p, str, bg_model
					);
				
				// Check if this score meets the threshold
				if(log_diff_score > thresh_score){

					// Add to the motif list
					motif_list.add(FindHighScoringRegionMotifs.get_new_FoundMotif(
						name, chr, p, str, log_diff_score, pwm[0].length
						));
					
				}
				
			}
		
		}
		
		// Sort the motif_list by motif score
		Collections.sort(motif_list);
		
		// Return the motif_list
		return(motif_list);
		
	}
	
	public static void write_motif_output( 
		BufferedWriter output, ArrayList<FoundMotif> motif_list 
		) throws IOException{
	
		// Iterate through the motif list
		for(FoundMotif fm : motif_list){
			
			// Write the output
			output.write(fm.getOutput() + "\n");
				
		}
			
	}
	
	

	public static void main(String[] args) throws IOException {

		// Set the input file name
		String feature_file_name = args[0]; 
					
		// Set the output file name
		String output_file_name = args[1];
		
		// Set the pwm file name
		String pwm_file_name = args[2];
			
		// Set the bg file name
		String genome_bg_file_name = args[3];

		// Set the log2 threshold score
		double thresh_score = Integer.parseInt(args[4]);

		// Set the genome fasta header
		String genome_fasta_header = args[5];
		
		// Set the feature_file_type
		boolean isSummit = Boolean.parseBoolean(args[6]);
		
		// Set the win if isPeak is false
		int win = 0;
		if(isSummit){
			win = Integer.parseInt(args[7]);
		}
		
		
		///////////////////////////////////////////////
		
		// Open the input buffer
		BufferedReader input = new BufferedReader(new FileReader(feature_file_name));

		// Read in the header
		String line = input.readLine();
	
		// Open the output buffer
		BufferedWriter output = new BufferedWriter(new FileWriter(output_file_name));
		
		// Write the header
		output.write("name,chr,pos,strand,motif_start,motif_end,motif_score,peak_start,peak_end\n");
		
		// Set the chr and chr_seq
		String chr = "";
		String chr_seq = "";

		// Read in the acs pwm
		double[][] tf_pwm = PWMFunctions.read_in_pwm(pwm_file_name);
			
		// Get the background model
		HashMap<String, Double> bg_model = PWMFunctions.read_in_bg_model(genome_bg_file_name);

		
		///////////////////////////////////////////////////////////////////////////
		
		// Initialize the feature_start and feature_end
		int feature_start;
		int feature_end;
				
		while((line = input.readLine()) != null){
			
			// Split the input
			String[] line_arr = line.split(",");
			
			// Get the parameters
			String name = line_arr[0];
			String feature_chr = line_arr[1];
			
			// If peak, get the feature_start and end from the feature_file
			if(isSummit){
				feature_start = Integer.parseInt(line_arr[2]) - win;
				feature_end = Integer.parseInt(line_arr[2]) + win;
			}else{
				feature_start = Integer.parseInt(line_arr[2]);
				feature_end = Integer.parseInt(line_arr[3]);	
			}
						
			// Read in the chr_seq if necessary
			if(!feature_chr.equals(chr)){
				chr = feature_chr;
				chr_seq =  GenomeFunctions.read_in_chr_fasta_file(feature_chr, genome_fasta_header);
			}
			
			// Set the storage ArrayList
			ArrayList<FoundMotif> motif_list =	get_motif_list(
				name, chr, feature_start, feature_end, 
				tf_pwm, chr_seq, bg_model, thresh_score
				);
			
			// If the motif_list is not empty, write the output
			if(!motif_list.isEmpty()){
			
				write_motif_output(output, motif_list);
			
			}
							
		}
		
		// Close the input buffer
		input.close();
		
		// Close the output buffer
		output.close();
		
	}

}