/*
* Copyright 2012-2015 Broad Institute, Inc.
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
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.utils.report;

import org.broadinstitute.gatk.utils.commandline.Gatherer;
import org.broadinstitute.gatk.utils.exceptions.UserException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.List;

public class GATKReportGatherer extends Gatherer {
    @Override
    public void gather(List<File> inputs, File output) {
        //Combines inputs GATKReport to one output

        PrintStream o;
        try {
            o = new PrintStream(output);
        } catch (FileNotFoundException e) {
            throw new UserException(String.format("File %s to be output by GATKReportGatherer function was not found", output));
        }

        GATKReport current = new GATKReport();
        boolean isFirst = true;
        for (File input : inputs) {
            if (isFirst) {
                current = new GATKReport(input);
                isFirst = false;
            } else {
                current.concat(new GATKReport(input));
            }
        }

        current.print(o);
        o.close();
    }
}
