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

package edu.berkeley.cs.amplab.avocado.input

import edu.berkeley.cs.amplab.adam.avro.ADAMRecord
import edu.berkeley.cs.amplab.avocado.stats.AvocadoConfigAndStats
import org.apache.commons.configuration.{HierarchicalConfiguration, SubnodeConfiguration}
import org.apache.spark.rdd.RDD

object Input {

  // all our input stages
  val stages = List(AlignedReadsInputStage,
                    SnapInputStage)

  /**
   * Builds the input stage that corresponds to the given stage name, and returns the read data
   * that the stage provides. The input stage name to use is collected from the provided
   * configuration.
   *
   * @param inputPath Path to input read data.
   * @param config Configuration file containing the necessary data.
   * @param stats Global stat and configuration data.
   * @return Returns an RDD of read data.
   */
  def apply (inputPath: String,
             config: HierarchicalConfiguration,
             stats: AvocadoConfigAndStats): RDD[ADAMRecord] = {
    // get input stage to use; if none is specified, default to input being aligned reads
    val stageName: String = config.getString("inputStage", "AlignedReads")

    val stage = stages.find(_.stageName == stageName)
    
    stage match {
      case Some(s: InputStage) => {
        val stageConfig = config.configurationAt(stageName)
        
        s.apply(inputPath, stageConfig, stats)
      }
      case None => {
        throw new IllegalArgumentException("No input stage with name: " + stageName)
      }
    }
  }
}
