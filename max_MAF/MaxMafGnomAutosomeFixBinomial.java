import java.awt.Frame;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

public class MaxMafGnomAutosomeFixBinomial {
    
    //ueses a continuous Poisson probability to determine if an allele is rare enough based on alt count and total count
    
// /////////////////////////////////////////////////////// Poisson Distribution Engine "frequencyChecker" ////////////////////////////////////////////////
    
    public static double calcBinomCoeff(double n, int k) {
    
        double binomCoeff = 1.0;
        for (int i = 1; i <= k; i++) {
            binomCoeff = binomCoeff * ((n+1-i)/i);
        }
        return binomCoeff;
    }
    public static boolean frequencyChecker(String alleleCount, String coverage, double threshold) { // method for testing if any one group freq > threshold
        int altCount = Integer.parseInt(alleleCount); 
        int totCount = Integer.parseInt(coverage);
        
        if (totCount > 700) {                   // a totCount>700 is solid, no need for Possion, just use simple division for true freq estimate
            double rate = ((double) altCount) / ((double) totCount);
            if(rate > 1.2*threshold) {            //is the rate greater than the threshold? make T/F call
                return true;
            }
            return false;
        }
        
        //prob is the cumulative probability that the actual value is less than or equal to the threshold
        double prob = 0.0; 
        for (int i = 0; i <= altCount; i++) {          // iterate up to count for this total count to be threshold
            prob += calcBinomCoeff(totCount,i)*Math.pow(threshold,i)*Math.pow((1.0-threshold),(totCount-i));
        }
        
        if ((1.0 - prob) <= 0.05) {      //is the probability that the actual value is less than the threshold low enough?
            return true;        //if the probability is below, yes reject this - the value is probably over the threshold with 95% confidence!
        }
        return false;           // 95% chance the true probability is OK and this is probably a rare allele 
    }

    // ///////////////////////////////////////////// main program //////////////////////////////////////////////////////////////////////////////////////////////
    public static void main(String[] args) throws IOException {
        JFileChooser browseTo = new JFileChooser();
        JOptionPane.showMessageDialog(new Frame("Input prompt"), "Please select an input GnomAD Master file to add MAF fields too.");               
        browseTo.showOpenDialog(new JPanel());
        File tempFile = browseTo.getSelectedFile();
        
        BufferedReader exacReader = new BufferedReader(new FileReader(tempFile));
        PrintWriter writer = new PrintWriter(tempFile.getParent() + "\\" + tempFile.getName().substring(0, tempFile.getName().indexOf(".")) + "_MAF.txt");

        //exacReader.readLine();                            // optional skips a line if ##headder is present. not needed if first line is headder
        String curLine = exacReader.readLine();             //read a new line in with pop data for that variant allele
        String[] curLineSplit = curLine.split("\t");            // splits the input line into separate variables
        List<String> headers = Arrays.asList(curLineSplit);     // get the headder lines for this array of ethnic groups
        
        ArrayList<Integer> temp = new ArrayList<Integer>();     // make an array list of integers of alt counts and chrom counts
        temp.add(headers.indexOf("AC_AFR"));    // create the African data alt allele count headder
        temp.add(headers.indexOf("AN_AFR"));    // create the African data total allele count (half anyway) headder
        
        temp.add(headers.indexOf("AC_AMR"));    // create the AFR alt allele count headder
        temp.add(headers.indexOf("AN_AMR"));    // dreate the AFR total allele count (half anyway) headder

        temp.add(headers.indexOf("AC_ASJ"));    // East Asian headder for alt allele count 
        temp.add(headers.indexOf("AN_ASJ"));    // East Asian headder for total allele count (half anyway)

        temp.add(headers.indexOf("AC_EAS"));    // South Asian headder for alt allele count 
        temp.add(headers.indexOf("AN_EAS"));    // South Asian headder for total allele count (half anyway)

        temp.add(headers.indexOf("AC_FIN"));    // Nonfin European headder for alt allele count 
        temp.add(headers.indexOf("AN_FIN"));    // Nonfin Euro headder for total allele count (half anyway)

        temp.add(headers.indexOf("AC_NFE"));    // Latino headder for alt allele count 
        temp.add(headers.indexOf("AN_NFE"));    // Latino headder for alt allele count (half anyway)
        
        temp.add(headers.indexOf("AC_SAS"));    // South Asian headder for alt allele count 
        temp.add(headers.indexOf("AN_SAS"));    // South Asian headder for alt allele count (half anyway)

        temp.add(headers.indexOf("AC"));        // Total headder for alt allele count 
        temp.add(headers.indexOf("AN"));        // Total headder for total allele count (half anyway)
        
        if(temp.contains(-1)) {                 // inappropriate print if something is wrong
            System.out.println("fuck");
            System.out.println(temp.indexOf(-1));
            System.gc();
            System.exit(0);
        }
        
        ArrayList<int[]> pops = new ArrayList<int[]>();     // load the array with the integers corresponding to the
        for (int i = 0; i < temp.size() / 2; i++) {     // correct sequence of variables from the ethnic data
            int[] tempArray = new int[2];               // these will be loaded into the inv poisson distribution
            tempArray[0] = temp.get(i*2);               // calculator method (frequencyChecker() for testing to
            tempArray[1] = temp.get(i*2 + 1);           // see if any have a freq > threshold at P<0.05 which
            pops.add(tempArray);                        // means the variant is too frequent to cause disease!
        }
        
        writer.println(curLine + "\t" +  "GnAC_Max_MAF>1%" + "\t" + "GnAC_Max_MAF>2%"); //writes out the headders for conclusions
        
        String f1;
        String f2;
        while (exacReader.ready()) {
            curLine = exacReader.readLine();
            curLineSplit = curLine.split("\t");
            f1 = "false";
            f2 = "false";
            
            for (int[] i: pops) {
                if (Integer.parseInt(curLineSplit[i[0]]) == 0){
                    // skip AC = 0 subpopulations
                    continue;
                }
                if (frequencyChecker(curLineSplit[i[0]], curLineSplit[i[1]], 0.01)) {
                    f1 = "true";
                }
                if (frequencyChecker(curLineSplit[i[0]], curLineSplit[i[1]], 0.02)) {
                    f2 = "true";
                }               
            }
            writer.println(curLine + "\t" +  f1 + "\t" + f2);
        }
        
        writer.close();
        exacReader.close();
        JOptionPane.showMessageDialog(new JPanel(), "Completed");
        System.gc();
        System.exit(0);
    }
}

