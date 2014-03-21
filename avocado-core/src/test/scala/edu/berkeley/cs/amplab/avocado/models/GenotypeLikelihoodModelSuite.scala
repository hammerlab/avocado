/*
 * Copyright (c) 2013. Regents of the University of California
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

package edu.berkeley.cs.amplab.avocado.models

import edu.berkeley.cs.amplab.adam.avro.{ADAMPileup, Base}
import org.scalatest.FunSuite
import scala.math.abs
import org.broad.tribble.gelitext.DiploidGenotype
import edu.berkeley.cs.amplab.adam.util.PhredUtils

class GenotypeLikelihoodModelSuite extends FunSuite {

  val floatingPointingThreshold = 1e-6

  def assertAlmostEqual(a: Double, b: Double, epsilon : Double = floatingPointingThreshold) {
    assert((a * 0.99 < b && a * 1.01 > b) ||
           abs(a - b) < epsilon)
  }

  test("score genotype for single sample, all bases ref") {
    val model = new DiploidGenotypeLikelihoodModel(2)

    val errorPhred30 = PhredUtils.phredToErrorProbability(30)
    val errorPhred40 = PhredUtils.phredToErrorProbability(40)

    val pl = List(ADAMPileup.newBuilder()
                    .setReadBase(Base.C)
                    .setReferenceBase(Base.C)
                    .setMapQuality(30)
                    .setSangerQuality(30)
                    .setCountAtPosition(1)
                    .build(),
                  ADAMPileup.newBuilder()
                    .setReadBase(Base.C)
                    .setReferenceBase(Base.C)
                    .setMapQuality(40)
                    .setSangerQuality(40)
                    .setCountAtPosition(1)
                    .build(),
                  ADAMPileup.newBuilder()
                    .setReadBase(Base.C)
                    .setReferenceBase(Base.C)
                    .setMapQuality(40)
                    .setSangerQuality(30)
                    .setCountAtPosition(1)
                    .build()
    )

    val hetLikelihood =  ((1 - errorPhred30) + errorPhred30 / 3) * ((1 - errorPhred40) + errorPhred40 / 3) * ((1 - errorPhred30) + errorPhred30 / 3) / 8.0
    val altLikelihood =  (2 * errorPhred30 / 3) * (2 * errorPhred40 / 3) * (2 * errorPhred30 / 3) / 8.0
    val expectedLikelihoods  =  scala.collection.mutable.Map.empty[DiploidGenotype, Double]

    expectedLikelihoods +=  DiploidGenotype.CC -> (2 * ((1 - errorPhred30) * 2 * (1 - errorPhred40) * 2 * (1 - errorPhred30))) / 8.0

    expectedLikelihoods += DiploidGenotype.AC -> hetLikelihood
    expectedLikelihoods += DiploidGenotype.CG -> hetLikelihood
    expectedLikelihoods += DiploidGenotype.CT -> hetLikelihood


    DiploidGenotype.values.filterNot(expectedLikelihoods.contains).map(gt => expectedLikelihoods.put(gt, altLikelihood ))

    val scored = model.computeLikelihoods(pl)
    scored.foreach( l => assertAlmostEqual(l._2, expectedLikelihoods(l._1)))
  }

//  test("score genotype for single sample, mix of ref/non-ref bases") {
//    val call = new DiploidGenotypeLikelihoodModel(2)
//
//    val pl = List(ADAMPileup.newBuilder()
//                    .setReadBase(Base.C)
//                    .setReferenceBase(Base.C)
//                    .setMapQuality(30)
//                    .setSangerQuality(30)
//                    .setCountAtPosition(1)
//                    .build(),
//                  ADAMPileup.newBuilder()
//                    .setReadBase(Base.C)
//                    .setReferenceBase(Base.C)
//                    .setMapQuality(40)
//                    .setSangerQuality(40)
//                    .setCountAtPosition(1)
//                    .build(),
//                  ADAMPileup.newBuilder()
//                    .setReadBase(Base.A)
//                    .setReferenceBase(Base.C)
//                    .setMapQuality(40)
//                    .setSangerQuality(30)
//                    .setCountAtPosition(1)
//                    .build()
//                  )
//    //TODO(arahuja) test is not correct as multiplying probabilities is not exactly the same as averaging qual scores
//    val expected = List((8.0 * ((0.999 * 0.999) * (0.9999 * 0.9999) * (1.0 - 0.999 * 0.9999))) / 8.0,
//
//      ((0.999 * 0.999 + (1.0 - 0.999 * 0.999)) *  (0.9999 * 0.9999 + (1.0 - 0.9999 * 0.9999))
//        *  ((1.0 - 0.999 * 0.9999) + (0.999 * 0.9999) )) / 8.0,
//      (8.0 * (1.0 - 0.999 * 0.999) * (1.0 - 0.9999 * 0.9999) * (0.999 * 0.9999)) / 8.0).reverse
//
//    val scored = call.scoreGenotypeLikelihoods(pl)
//
//    for (i <- 0 to 2) {
//      assertAlmostEqual(expected(i), scored(i), 1e-3)
//    }
//  }
//
//  test("score genotype for single sample, all bases non-ref") {
//    val call = new PileupCallSimpleSNP(2)
//
//    val pl = List(ADAMPileup.newBuilder()
//                    .setReadBase(Base.A)
//                    .setReferenceBase(Base.C)
//                    .setMapQuality(30)
//                    .setSangerQuality(30)
//                    .setCountAtPosition(1)
//                    .build(),
//                  ADAMPileup.newBuilder()
//                    .setReadBase(Base.A)
//                    .setReferenceBase(Base.C)
//                    .setMapQuality(40)
//                    .setSangerQuality(40)
//                    .setCountAtPosition(1)
//                    .build(),
//                  ADAMPileup.newBuilder()
//                    .setReadBase(Base.A)
//                    .setReferenceBase(Base.C)
//                    .setMapQuality(40)
//                    .setSangerQuality(30)
//                    .setCountAtPosition(1)
//                    .build()
//                  )
//
//    val expected = List(8.0 * (1.0 - 0.999 * 0.999) * (1.0 - 0.9999 * 0.9999) * (1.0 - 0.999 * 0.9999) / 8.0,
//                        ((0.999 * 0.999 * 0.9999 * 0.9999 * 0.999 * 0.9999) +
//                         (1.0 - 0.999 * 0.999) * (1.0 - 0.9999 * 0.9999) * (1.0 - 0.999 * 0.9999)) / 8.0,
//                        (8.0 * (0.999 * 0.999 * 0.9999 * 0.9999 * 0.999 * 0.9999)) / 8.0).reverse
//
//
//    val scored = call.scoreGenotypeLikelihoods(pl)
//
//    for (i <- 0 to 2) {
//      assertAlmostEqual(expected(i), scored(i))
//    }
//  }

}
