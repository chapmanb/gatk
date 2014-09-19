/**
   From Daniel Klevebring's GATK Repository: https://github.com/dakl/gatk
*/
package org.broadinstitute.gatk.engine.walkers.dakl.utils.MultipleTestingCorrection;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

/**
 * Created with IntelliJ IDEA.
 * User: dankle
 * Date: 2014-02-18
 * Time: 20:51
 * To change this template use File | Settings | File Templates.
 */
public class MultipleTestingCorrection {


	/*--------------------------------------------------------------
      FIELDS.
      --------------------------------------------------------------*/
    /** significance level.*/
    private static String alpha;
    /** hashmap with the test results as values and the GO labels as keys.*/
    private static HashMap map;
    /** type of correction*/
    private static String type ;
    /** hashmap with the results (adjusted p-values) as values and the GO labels as keys.*/
    private static HashMap<String, Double> correctionMap;

	/*--------------------------------------------------------------
      CONSTRUCTOR.
      --------------------------------------------------------------*/

    public MultipleTestingCorrection (String alpha, HashMap map){
        this.alpha = alpha ;
        this.map = map ;
        this.correctionMap = null ;
    }



	/*--------------------------------------------------------------
		METHODS.
      --------------------------------------------------------------*/

    /**
     * method that redirects the calculation of the multiple testing correction.
     */
    public void calculate(){

        HashSet goLabelsSet = new HashSet(map.keySet());

        Iterator iteratorGoLabelsSet = goLabelsSet.iterator();
        String [] pvalues = new String [map.size()];
        String [] goLabels = new String [map.size()];
        for(int i = 0; iteratorGoLabelsSet.hasNext(); i++){
            goLabels[i] = iteratorGoLabelsSet.next().toString();
            pvalues[i] = map.get(new String(goLabels[i])).toString();
        }

        String [] adjustedPvalues ;
        String [] sortedGOLabels ;

        BenjaminiHochbergFDR fdr = new BenjaminiHochbergFDR(pvalues, goLabels, alpha);
        fdr.calculate();
        adjustedPvalues = fdr.getAdjustedPvalues();
        sortedGOLabels = fdr.getOrdenedGOLabels();
        correctionMap = new HashMap<String,Double>();
        for (int i = 0; i < adjustedPvalues.length && i < sortedGOLabels.length; i++){
            correctionMap.put(new String(sortedGOLabels[i]), new Double(adjustedPvalues[i]));
        }

    }


	/*--------------------------------------------------------------
		GETTER.
      --------------------------------------------------------------*/

    /**
     * getter for the map of corrected p-values.
     *
     * @return HashMap correctionMap.
     */
    public HashMap getCorrectionMap(){
        return correctionMap;
    }

}
