package motif_processing;

import java.io.IOException;
import java.util.Arrays;

public class CMDLineMain {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
				
		// Get the program to execute
		String java_executable_program = args[0];
		
		// Get the command line arguments
		String[] cmd_line_args = Arrays.copyOfRange(args, 1, args.length);
				
		// Choose the program to run based on the args[0]
		if(java_executable_program.equals("FindHighScoringRegionMotifs")){
			System.out.println("Executing FindHighScoringRegionMotifs!");
			FindHighScoringRegionMotifs.main(cmd_line_args);
		}else if(java_executable_program.equals("FindMotifInRegion")){
			System.out.println("Executing FindMotifInRegion!");
			FindMotifInRegion.main(cmd_line_args);
		}else{
			System.out.println("Did not specify either 'FindHighScoringRegionMotifs' or " +
							   "'FindMotifInRegion' as first argument!");
		}

	}

}
