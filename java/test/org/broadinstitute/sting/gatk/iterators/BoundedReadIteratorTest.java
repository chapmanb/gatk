package org.broadinstitute.sting.gatk.iterators;

import static junit.framework.Assert.fail;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import net.sf.picard.reference.ReferenceSequenceFile;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.datasources.shards.Shard;
import org.broadinstitute.sting.gatk.datasources.shards.ShardStrategy;
import org.broadinstitute.sting.gatk.datasources.shards.ShardStrategyFactory;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMDataSource;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SimpleDataSourceLoadException;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.IndexDrivenSAMDataSource;
import org.broadinstitute.sting.gatk.Reads;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;



/*
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * @author aaron
 * @version 1.0
 * @date Apr 14, 2009
 * <p/>
 * Class BoundedReadIteratorTest
 * <p/>
 * tests for the bounded read iterator.
 */
public class BoundedReadIteratorTest extends BaseTest {

    /** the file list and the fasta sequence */
    private List<File> fl;
    private ReferenceSequenceFile seq;

    /**
     * This function does the setup of our parser, before each method call.
     * <p/>
     * Called before every test case method.
     */
    @Before
    public void doForEachTest() throws FileNotFoundException {
        fl = new ArrayList<File>();

        // sequence
        seq = new IndexedFastaSequenceFile(new File(seqLocation + "/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta"));
        GenomeLocParser.setupRefContigOrdering(seq.getSequenceDictionary());
    }


    /** Test out that we can shard the file and iterate over every read */
    @Test
    public void testBounding() {
        logger.warn("Executing testBounding");
        // total reads expected
        final int expected = 20;
        // bound by ten reads
        BoundedReadIterator iter = new BoundedReadIterator(new testIterator(), expected);

        int count = 0;
        for (SAMRecord rec: iter) {
            count++;
        }

        Assert.assertEquals(expected,count);
    }
}

class testIterator implements StingSAMIterator {
    SAMFileHeader header;
    testIterator() {
        header = ArtificialSAMUtils.createArtificialSamHeader(1,1,2000);
    }
    /**
     * Gets source information for the reads.  Contains information about the original reads
     * files, plus information about downsampling, etc.
     *
     * @return
     */
    public Reads getSourceInfo() {
        return null;
    }

    public void close() {

    }

    public boolean hasNext() {
        return true;
    }

    public SAMRecord next() {
        return ArtificialSAMUtils.createArtificialRead(header,"blah",0,1,100);
    }

    public void remove() {
    }

    public Iterator<SAMRecord> iterator() {
        return this;
    }
}