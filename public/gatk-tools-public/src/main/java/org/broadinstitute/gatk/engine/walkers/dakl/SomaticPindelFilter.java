/**
   From Daniel Klevebring's GATK Repository: https://github.com/dakl/gatk
*/

package org.broadinstitute.gatk.engine.walkers.dakl;

import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.gatk.engine.walkers.dakl.utils.FisherExact;
import org.broadinstitute.gatk.engine.walkers.dakl.utils.MultipleTestingCorrection.*;

import java.util.*;

/**
 *
 * Filters pindel output based on coverage and Fisher's exact test (BH corrected).
 *
 * 1. Quality filters
 *   Tumor read depth > 25 and < 1000 (changeable from cli)
 *   Normal read depth > 25 and < 1000 (changeable from cli)
 *   Microhomology around indel no longer than 2 bases (ie HOMLEN < abs(SVLEN) + 2)
 *
 * 2. Test for difference between tumor allele fraction and normal
 *   Perform a one-tailed Fisher's exact test on each variant
 *   Adjust for multiple testing using Benjamini Hochberg
 *
 * Variants are written in VCF format, where the INFO field PV is added (with the adjusted p value) and
 * the FORMAT field FA with the allele fractions
 *
 * Created with IntelliJ IDEA.
 * User: dankle
 * Date: 2014-01-31
 * Time: 20:05
 * To change this template use File | Settings | File Templates.
 */
public class SomaticPindelFilter extends RodWalker<Integer,Integer> {
    @Output(doc="File to which variants should be written")
    protected VariantContextWriter vcfWriter = null;

    @Input(fullName="variant", shortName = "V", doc="Select variants from this VCF file",
            required=true)
    public RodBinding<VariantContext> variants;

    @Input(fullName="tumorid", shortName = "TID", doc="Sample Name of the tumor", required = true)
    public String tumorId;

    @Input(fullName="normalid", shortName = "NID", doc="Sample Name of the normal", required = true)
    public String normalId;

    @Input(fullName="minCoverageNormal", shortName = "MIN_DP_N", doc="Minimum depth required in the normal", required = false)
    public int MIN_DP_N = 25;

    @Input(fullName="minCoverageTumor", shortName = "MIN_DP_T", doc="Minimum depth required in the tumor", required = false)
    public int MIN_DP_T = 25;

    @Input(fullName="maxCoverageNormal", shortName = "MAX_DP_N", doc="Maximum depth required in the normal", required = false)
    public int MAX_DP_N = 1000;

    @Input(fullName="maxCoverageTumor", shortName = "MAX_DP_T", doc="Maximum depth required in the tumor", required = false)
    public int MAX_DP_T = 1000;

    @Input(fullName="maxNormalFraction", shortName = "MAX_N_FRAC", doc="Maximum fraction allowed in the normal", required = false)
    public double MAX_N_FRAC = .15;

    @Input(fullName="AdjustedPValueCutoff", shortName = "ADJ_P_CUTOFF", doc="Cutoff for the adjusted p values", required = false)
    public double ADJ_P_CUTOFF = 0.05;

    private FisherExact fisher = new FisherExact(100000);
    private ArrayList<VariantContext> candidateVariants;

    /**
     * Set up the VCF writer, the sample expressions and regexs, and the JEXL matcher
     */
    public void initialize() {
        // Set up new header including PV INFO field
        Map<String, VCFHeader> vcfRods = GATKVCFUtils.getVCFHeadersFromRods(getToolkit());
        Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), true);
        headerLines.add(new VCFHeaderLine(VCFHeader.SOURCE_KEY, "SomaticPindelFilter"));
        TreeSet<String> vcfSamples = new TreeSet<String>(SampleUtils.getSampleList(vcfRods, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE));
        VCFInfoHeaderLine vcfInfoHeaderLine = new VCFInfoHeaderLine("PVP", 1, VCFHeaderLineType.Float, "P value after adjustment for multiple testing using Benjamini Hochberg's method for controlling FDR (from SomaticPindelFilter).");
        VCFFormatHeaderLine vcfFormatHeaderLine = new VCFFormatHeaderLine("FA", 1, VCFHeaderLineType.Float, "Allele fraction of the alternate allele with regard to reference");
        headerLines.add( vcfInfoHeaderLine );
        headerLines.add( vcfFormatHeaderLine );
        VCFHeader vcfHeader = new VCFHeader(headerLines, vcfSamples);
        vcfWriter.writeHeader(vcfHeader);

        // set up list of candidate variants. These variants will be tested with Fishers exact test to check for AF T and N
        candidateVariants = new ArrayList<VariantContext>();
    }

        @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        Collection<VariantContext> VCs = tracker.getValues(variants, context.getLocation());
        for ( VariantContext vc : VCs ){
            /*
            * Keep only those variants that have HOMLEN < abs(SVLEN+3)
            * */
            int homLen = vc.getAttributeAsInt("HOMLEN", 100);
            int svLen  = vc.getAttributeAsInt("SVLEN", 0);
            int t_dp   = vc.getGenotype(tumorId).getAD()[0] + vc.getGenotype(tumorId).getAD()[1];
            int n_dp   = vc.getGenotype(normalId).getAD()[0] + vc.getGenotype(normalId).getAD()[1];
            double n_frac = (double)vc.getGenotype(normalId).getAD()[1] / (double)n_dp;
            int absoluteSvLen = (svLen < 0) ? -svLen : svLen;
            // keep variants if no microhomology is present AND
            // coverage is in correct range
            // and normal fraction is <= MAX_N_FRAC
            if( homLen < absoluteSvLen + 2 && t_dp > MIN_DP_T && t_dp < MAX_DP_T && n_dp > MIN_DP_N && n_dp < MAX_DP_N && n_frac <= MAX_N_FRAC){
                candidateVariants.add(vc);
            }
        }
        return 0;
    }

    @Override
    public Integer reduceInit() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public Integer reduce(Integer value, Integer sum) {
        return value + sum;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public void onTraversalDone(Integer result) {
        HashMap<String,Double> pvalues = new HashMap<String,Double>();
        for( VariantContext vc : candidateVariants ){
            /*    2x2 tab for fisher's exact test
            *             REF    ALT
            *      TUMOR    a      b
            *      NORMAL   c      d
            *
            * */
            int a = vc.getGenotype(tumorId).getAD()[0];
            int b = vc.getGenotype(tumorId).getAD()[1];
            int c = vc.getGenotype(normalId).getAD()[0];
            int d = vc.getGenotype(normalId).getAD()[1];
            double p = fisher.getLeftTailedP(a,b,c,d);
            pvalues.put( vc.getChr()+":"+vc.getStart()+":"+vc.getReference()+":"+vc.getAltAlleleWithHighestAlleleCount(), p);
        }
        MultipleTestingCorrection mtc = new MultipleTestingCorrection(Double.toString(ADJ_P_CUTOFF), pvalues);
        mtc.calculate();

        for( VariantContext vc : candidateVariants ){
            double padj = Double.parseDouble(mtc.getCorrectionMap().get( vc.getChr()+":"+vc.getStart()+":"+vc.getReference()+":"+vc.getAltAlleleWithHighestAlleleCount() ).toString() );

            if( padj < ADJ_P_CUTOFF ){ // write the variant if it's significant
                Map attributes = vc.getAttributes();
                GenotypesContext genotypesContext = vc.getGenotypes();
                ArrayList<Genotype> newGenotypes = new ArrayList<Genotype>();
                for( Genotype g : genotypesContext ){
                    double allele_frac = (double)g.getAD()[1]/((double)g.getAD()[0] + (double)g.getAD()[1]);
                    Genotype newGenotype = new GenotypeBuilder(g).attribute("FA", allele_frac).make();
                    newGenotypes.add(newGenotype);
                }
                VariantContext newVariantContext = new VariantContextBuilder(vc).genotypes(newGenotypes).attribute("PVP", padj).make();
                vcfWriter.add(newVariantContext);
            }
        }
    }

}

