package edu.berkeley.cs.amplab.avocado.models

import edu.berkeley.cs.amplab.adam.avro.ADAMPileup
import org.apache.spark.rdd.PairRDDFunctions
import edu.berkeley.cs.amplab.adam.util.PhredUtils
import edu.berkeley.cs.amplab.adam.util.Base
import org.broad.tribble.gelitext.DiploidGenotype

abstract class GenotypeLikelihoodModel(val ploidy : Int = 2) {

  object Genotype extends Enumeration {
    type Genotype = Value
  }

  val BASES = List('A','C', 'T', 'G')
  val N_BASES = BASES.size
  val PCR_ERROR = math.log10(10e-4)

  protected def computeBaseGenotypeLikelihood( base : ADAMPileup, genotype : Seq[Char]) : Double = {
    genotype.map(refBase => computeBaseLikelihood(base, refBase)).sum
  }

  protected def computeBaseGenotypeLogLikelihood( base : ADAMPileup, genotype : Seq[Char]) : Double = {
    math.log10(genotype.map(refBase => computeBaseLikelihood(base, refBase)).sum)
  }

  protected def computeBaseLikelihood( base : ADAMPileup, refBase : Char) : Double = {
    val error = PhredUtils.phredToErrorProbability(base.getSangerQuality)
    return if (base.getReadBase.toString == refBase.toString) (1 - error) else error / (N_BASES - 1)
  }

  protected def computeBaseLogLikelihood( base : ADAMPileup, refBase : Char) : Double = {
    val error = PhredUtils.phredToErrorProbability(base.getSangerQuality)
    return if (base.getReadBase.toString == refBase.toString) math.log10(1 - error) else math.log10(error)
  }

  private def computePairedBaseLikelihood( firstRead : ADAMPileup , secondRead : ADAMPileup, refBase : Char) : Double = {

    def fragmentLikelihood(firstRead : ADAMPileup , secondRead : ADAMPileup, refBase : Char, fragmentBase : Char) : Double = {
      val baseLikelihood = computeBaseLikelihood(firstRead, fragmentBase) * computeBaseLikelihood(secondRead, fragmentBase)
      if (refBase == fragmentBase) baseLikelihood * (1 - PCR_ERROR) else baseLikelihood * PCR_ERROR/ (N_BASES - 1)
    }

    BASES.map(fragmentLikelihood(firstRead, secondRead, refBase, _)).sum
  }

  def computeLikelihoods(pileup: List[ADAMPileup]): Map[DiploidGenotype, Double]

}

class DiploidGenotypeLikelihoodModel(override val ploidy : Int = 2) extends GenotypeLikelihoodModel(ploidy) {

  def computeLikelihoods(pileup: List[ADAMPileup]): Map[DiploidGenotype, Double] = {

    def computeGenotypeLikelihood(pileup : List[ADAMPileup], gt : DiploidGenotype) : Double = {
      pileup.map( p => computeBaseGenotypeLikelihood(p, gt.toString)).reduce( _ * _ )
    }

    val k = pileup.map(_.getCountAtPosition).reduce(_ + _)
    val denom = math.pow(ploidy, k.toDouble)
    DiploidGenotype.values.map(gt => ( gt , (computeGenotypeLikelihood(pileup, gt) / denom))).toMap
  }

  def computeGenotypeLogLikelihood(pileup : List[ADAMPileup], gt : DiploidGenotype) : Double = {
    pileup.map( p => computeBaseGenotypeLogLikelihood(p, gt.toString)).reduce( _ + _ )
  }

}