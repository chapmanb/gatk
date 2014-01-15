/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.queue.extensions.gatk

import org.broadinstitute.sting.queue.function.RetryMemoryLimit
import org.broadinstitute.sting.queue.function.scattergather.GatherFunction


/**
 *
 * Currently this is the default gather for VCFs.
 * One can set a specific gatherer to use by adding @Gather before any output argument.
 * For example (used to be part of UG):
 *           @Gather(className = "org.broadinstitute.sting.queue.extensions.gatk.CatVariantsGatherer")
 *           @Output(doc="File to which variants should be written",required=true)
 *           protected VariantContextWriter writer = null;
 */
class CatVariantsGatherer extends CatVariants with GatherFunction with RetryMemoryLimit{
  this.assumeSorted = true

  private lazy val originalGATK = this.originalFunction.asInstanceOf[CommandLineGATK]

  override def freezeFieldValues() {
    this.reference = originalGATK.reference_sequence
    this.variant = this.gatherParts.zipWithIndex map { case (input, index) => new TaggedFile(input, "input"+index) }
    this.outputFile = this.originalOutput
    this.assumeSorted = true
    this.variant_index_type = originalGATK.variant_index_type
    this.variant_index_parameter = originalGATK.variant_index_parameter

    super.freezeFieldValues()
  }



}
