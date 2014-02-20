/*
 * Copyright (c) 2014. Regents of the University of California
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package edu.berkeley.cs.amplab.avocado.stats

import org.apache.spark.SparkContext
import org.apache.spark.rdd.RDD
import edu.berkeley.cs.amplab.adam.avro.{ADAMRecord, ADAMNucleotideContig}
import edu.berkeley.cs.amplab.adam.models.SequenceDictionary
import edu.berkeley.cs.amplab.adam.rdd.AdamContext._
import edu.berkeley.cs.amplab.adam.rdd.AdamRDDFunctions

class AvocadoConfigAndStats (val sc: SparkContext,
                             val debug: Boolean, 
                             reference: RDD[ADAMNucleotideContig]) {
  var inputDataset: RDD[ADAMRecord] = null
  
  lazy val coverage = ScoreCoverage(inputDataset)

  lazy val contigLengths = GetReferenceContigLengths(reference)
  
  lazy val referenceSeq = reference.collect()

  lazy val sequenceDict = reference.adamGetSequenceDictionary()

  lazy val samplesInDataset = inputDataset.map(_.getRecordGroupSample)
    .distinct()
    .collect()

  /**
   * Attaches reads to the statistics bin.
   *
   * @note Reads must be attached before any of the read based statistics can be accessed.
   * If reads are not attached, a NullPointerException will be thrown.
   *
   * @param reads Reads to attach to config.
   *
   * @throws AssertionError Throws an assertion error if reads are attached multiple times.
   */
  def attachReads (reads: RDD[ADAMRecord]) {
    assert(inputDataset == null, "Cannot attach reads multiple times.")
    inputDataset = reads
  }
}
