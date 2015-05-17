package motif_processing;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Set;

import jbfunctions.GenomeFunctions;
import jbfunctions.PWMFunctions;


public class FindHighScoringRegionMotifs {


	
	public static FoundMotif get_new_FoundMotif(
		String tf_name, String chr, int pos, char str, double pwm_score, int pwm_width
		){
		
		// Set the motif_pos_start and motif_pos_end variables
		int motif_pos_start = pos;
		int motif_pos_end = motif_pos_start + pwm_width - 1;
		
		// Update the coordinates if the pwm is on the reverse strand
		if(str == '-'){
			motif_pos_end = pos;
			motif_pos_start = pos - pwm_width + 1;
		}
		
		// Return the FoundMotif
		return(new FoundMotif(tf_name, chr, motif_pos_start, motif_pos_end, str, pwm_score));
		
	}
	
	public static void write_motif_output(ArrayList<FoundMotif> motif_list,
										  String output_file_name
										 ) throws IOException{
		
		// Open the output buffer
		BufferedWriter output = new BufferedWriter(new FileWriter(output_file_name));
		
		// Write the header
		output.write("tf,chr,start,strand,motif_start,motif_end,pwm_score\n");
		
		// Iterate through each motif
		for(int i = 0; i < motif_list.size(); i++)
			output.write(motif_list.get(i).getOutput() + "\n");
		
		// Close the output buffer
		output.close();
		
	}
	
	public static HashMap<String, double[][]> read_in_tf_pwm_list(String dirPath) throws IOException{
		
		// Get the directory contents
		File dir = new File(dirPath);
		
		// Get the listOfFiles
		File[] listOfFiles = dir.listFiles(new FilenameFilter(){
			
			public boolean accept(File dir, String name){
				return(name.endsWith(".pwm"));
			}
		
		});	
		
		// Set up the HashMap
		HashMap<String, double[][]> pwm_hash = new HashMap<String, double[][]>();
		
		// Iterate through each file
		for(File f : listOfFiles){
			
			// Get the PWM String
			String pwm_string = f.getName();
			pwm_string = pwm_string.substring(0, pwm_string.length() - 4);
			
			// Get the PWM
			double[][] pwm_arr = PWMFunctions.read_in_pwm(f.getAbsolutePath());
			
			// Enter into the HashMap
			pwm_hash.put(pwm_string, pwm_arr);
			
		}
		
		// Return the HashMap
		return(pwm_hash);
	
	}
			
	public static void main(String[] args) throws IOException {
		
		// Read in the genomic location
		String chr = args[0];
		int start = Integer.parseInt(args[1]);
		int end = Integer.parseInt(args[2]);
		
		// Set the pwm_directory
		String pwm_dir = args[3];
		
		// Set the background model file
		String bg_model_file = args[4];
		
		// Set the output file
		String output_file_name = args[5];
		
		// Set the genome fasta directory
		String genome_fasta_dir = args[6];
				
		// Get the log2 threshold score
		double thresh_score = Integer.parseInt(args[7]);
				
		// Get the TF List
		HashMap<String, double[][]> pwm_hash = read_in_tf_pwm_list(pwm_dir);
		
		// Get the background model
		HashMap<String, Double> bg_model = PWMFunctions.read_in_bg_model(bg_model_file);
		
		// Get the sequence over the chromosome
		String chr_seq = GenomeFunctions.read_in_chr_fasta_file(chr, genome_fasta_dir);
		
		// Get the str array
		char[] str_arr = {'+', '-'};
		
		// Create the storage ArrayList for found motifs
		ArrayList<FoundMotif> motif_list = new ArrayList<FoundMotif>();
	
		// Get the list of TFs
		Set<String> tf_list = pwm_hash.keySet();
				
		// Iterate through each tf
		for(String tf : tf_list){
			
			// Read in the tf pwm
			double[][] tf_pwm = pwm_hash.get(tf);

			// Get the length of the pwm
			int pwm_length = tf_pwm[0].length;
							
			// Iterate over each strand
			for(char str : str_arr){
									
				for(int p = start; p <= end; p++){

					// Calculate the log-odds score
					double log_diff_score =	PWMFunctions.get_log_score(
												chr_seq, tf_pwm, p, str, bg_model
							   					);
					
					// Check if this score meets the threshold
					if(log_diff_score > thresh_score){

						// Add to the motif list
						motif_list.add(get_new_FoundMotif(tf, chr, p, str, log_diff_score, pwm_length));						
						
					}
				
				}
				
			}				
			
			// Sort the motifs by the PWM score
			Collections.sort(motif_list);
			
			// Write the output file
			write_motif_output(motif_list, output_file_name);

		}
		
	}
	
}
