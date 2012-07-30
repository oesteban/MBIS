/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkWeightedCentroidKdTreeImageGenerator_h
#define __itkWeightedCentroidKdTreeImageGenerator_h

#include <vector>

#include "itkKdTree.h"
#include "itkStatisticsAlgorithm.h"
#include "Numerics/Statistics/include/MaskedImageToListSampleAdaptor.h"

namespace itk
{
namespace Statistics
{
/** \class WeightedCentroidKdTreeImageGenerator
 *  \brief This class generates a KdTree object with centroid information.
 *
 * The KdTree object stores measurment vectors in a k-d tree structure
 * that is a binary tree. The partition value is the median value of one
 * of the k dimension (partition dimension). The partition dimension is
 * determined by the spread of measurement values in each dimension. The
 * partition dimension is the dimension has the widest spread. Our
 * implementation of k-d tree doesn't have any construction or insertion
 * logic. Users should use this class or the KdTreeGenerator class.
 *
 * This class is derived from the KdTreeGenerator class. The only
 * difference between this class and the KdTreeGenerator class is that
 * the nonterminal node type of this class is
 * KdTreeWeightedCentroidNonterminalNode and that of the
 * KdTreeGenerator is KdTreeNonterminalNode. Therefore, the public
 * interface is identical to each other. The nonterminal node generation
 * routines differ.
 *
 * To run this generator, users should provides the bucket size
 * (SetBucketSize method) and the input sample (SetSample method). The
 * Update method will run this generator. To get the resulting KdTree
 * object, call the GetOutput method.
 *
 * <b>Recent API changes:</b>
 * The static const macro to get the length of a measurement vector,
 * 'MeasurementVectorSize'  has been removed to allow the length of a measurement
 * vector to be specified at run time. It is now obtained from the sample set
 * as input. You may query this length using the function GetMeasurementVectorSize().
 *
 * \sa KdTree, KdTreeNode, KdTreeWeightedCentroidNonterminalNode,
 * KdTreeTerminalNode, KdTreeGenerator
 * \ingroup ITKStatistics
 */

template< class TImage, class TMask >
class ITK_EXPORT WeightedCentroidKdTreeImageGenerator: public Object
{
public:
  /** Standard class typedefs */
  typedef WeightedCentroidKdTreeImageGenerator  Self;
  typedef Object                                Superclass;
  typedef SmartPointer< Self >                  Pointer;
  typedef SmartPointer< const Self >            ConstPointer;

  /** Run-time type information (and related methods) */
  itkTypeMacro(WeightedCentroidKdTreeImageGenerator, Object);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  typedef TImage                                                          InputImageType;
  typedef TMask                                                           MaskImageType;
  typedef MaskedImageToListSampleAdaptor< TImage, TMask  >                SubsampleType;
  typedef typename SubsampleType::Pointer                                    SubsamplePointer;

  /** typedef alias for the source data container */
  typedef typename SubsampleType::MeasurementVectorType MeasurementVectorType;
  typedef typename SubsampleType::MeasurementType       MeasurementType;

  /** Typedef for the length of each measurement vector */
  typedef size_t MeasurementVectorSizeType;

  /** Typedef for the k-d tree */
  typedef KdTree< SubsampleType > KdTreeType;

  /** Type alias for the k-d tree type */
  typedef KdTreeType OutputType;

  /** Typedef for the smart pointer to the k-d tree */
  typedef typename KdTreeType::Pointer OutputPointer;

  /** Typedef for the k-d tree node type */
  typedef typename KdTreeType::KdTreeNodeType KdTreeNodeType;

  typedef KdTreeWeightedCentroidNonterminalNode< SubsampleType > KdTreeNonterminalNodeType;

  itkSetConstObjectMacro( InputImage, InputImageType);
  itkGetConstObjectMacro( InputImage, InputImageType);
  itkSetConstObjectMacro(  MaskImage, MaskImageType);
  itkGetConstObjectMacro(  MaskImage, MaskImageType);

  void SetBucketSize(size_t size);

  /** Returns the pointer to the generated k-d tree. */
  OutputPointer GetOutput()
  {
    return m_Tree;
  }

  /** Runs this k-d tree construction algorithm. */
  void Update();

  /** Get macro to get the length of the measurement vectors that are being
   * held in the 'sample' that is passed to this class */
  itkGetConstMacro(MeasurementVectorSize, size_t);

protected:
  /** Constructor */
  WeightedCentroidKdTreeImageGenerator();

  /** Destructor */
  virtual ~WeightedCentroidKdTreeImageGenerator() {}

  void PrintSelf(std::ostream & os, Indent indent) const;

  /** Nonterminal node generation routine */
  virtual KdTreeNodeType * GenerateNonterminalNode(size_t beginIndex,
                                                   size_t endIndex,
                                                   MeasurementVectorType
                                                   & lowerBound,
                                                   MeasurementVectorType
                                                   & upperBound,
                                                   size_t level);

  virtual KdTreeNodeType * GenerateTreeLoop(size_t beginIndex,
                                                   size_t endIndex,
                                                   MeasurementVectorType & lowerBound,
                                                   MeasurementVectorType & upperBound,
                                                   size_t level);

  SubsamplePointer GetSubsample() {  return m_Subsample;  }


private:
  WeightedCentroidKdTreeImageGenerator(const Self &); //purposely not implemented
  void operator=(const Self &);                  //purposely not implemented


  typename InputImageType::ConstPointer m_InputImage;
  typename MaskImageType::ConstPointer  m_MaskImage;

  /** Smart pointer to the internal Subsample object. This class needs
   * a Subsample object because the partitioning process involves sorting
   * and selection. */
  SubsamplePointer m_Subsample;

  /** The number of measurement vectors that can be stored in a terminal
   * node. */
  size_t m_BucketSize;

  /** Pointer to the resulting k-d tree. */
  OutputPointer m_Tree;

  /** Temporary lower bound for the TreeGenerationLoop */
  MeasurementVectorType m_TempLowerBound;

  /** Temporary upper bound for the TreeGenerationLoop */
  MeasurementVectorType m_TempUpperBound;

  /** Temporary mean for the TreeGenerationLoop */
  MeasurementVectorType m_TempMean;

  /** Length of a measurement vector */
  MeasurementVectorSizeType m_MeasurementVectorSize;

};  // end of class
} // end of namespace Statistics
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "WeightedCentroidKdTreeImageGenerator.hxx"
#endif

#endif
