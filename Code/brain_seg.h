// --------------------------------------------------------------------------------------
// File:          brain_seg.h
// Date:          Dec 14, 2011
// Author:        code@oscaresteban.es (Oscar Esteban)
// Version:       1.0 beta
// License:       GPLv3 - 29 June 2007
// Short Summary:
// --------------------------------------------------------------------------------------
//
// Copyright (c) 2012, code@oscaresteban.es (Oscar Esteban)
// with Signal Processing Lab 5, EPFL (LTS5-EPFL)
// and Biomedical Image Technology, UPM (BIT-UPM)
// All rights reserved.
//
// This file is part of MBIS
//
// MBIS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MBIS is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MBIS.  If not, see <http://www.gnu.org/licenses/>.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#ifndef BRAIN_SEG_H_
#define BRAIN_SEG_H_


#include <ostream>
#include <exception>
#include <cstdlib>
#include <fstream>
#include <iomanip>

#include <itkIntensityWindowingImageFilter.h>
#include <itkStatisticsImageFilter.h>
#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkThresholdImageFilter.h>
#include <itkLabelImageToLabelMapFilter.h>
#include <itkLabelMapMaskImageFilter.h>

#include <itkSubtractImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkMaskImageFilter.h>

#include <boost/unordered_map.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <boost/any.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>


#include "MultichannelParameterEstimator.h"
#include "MultispectralEMFilter.h"
#include "MultispectralKMeansFilter.h"
#include "MultispectralMLClassifierFilter.h"
#include "MultispectralMRFFilter.h"

using namespace std;
namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;

#define MRFSegmentation_VERSION_MAJOR 1 //@MRFSegmentation_VERSION_MAJOR@
#define MRFSegmentation_VERSION_MINOR 0 //@MRFSegmentation_VERSION_MINOR@

const unsigned int DEFAULT_NCLASSES         = 3;
const float DEFAULT_ATLAS_TH = 0.50;
const unsigned int DEFAULT_MAX_ITER         = 10;
const float DEFAULT_LAMBDA                  = 0.6;
const unsigned int DEFAULT_MRF_ITERATIONS   = 1;
const unsigned short NORM_MAX_INTENSITY     = 1000.0;
const std::string DEFAULT_MRF_ALG           = "expansion";

typedef float                                                                           PixelType;
typedef itk::Image< PixelType, 3 >                                                      InputComponentType;
typedef InputComponentType::Pointer                                                     InputComponentPointer;
typedef InputComponentType::ConstPointer                                                InputComponentConstPointer;

typedef float																			ChannelPixelType;
typedef itk::Image< ChannelPixelType, 3 >                                               ChannelImageType;
typedef ChannelImageType::Pointer                                                            ChannelPointer;
typedef ChannelImageType::ConstPointer                                                       ChannelConstPointer;
typedef itk::ImageFileReader< ChannelImageType >                                             ChannelReader;
typedef itk::CastImageFilter< ChannelImageType, InputComponentType >                         ChannelToComponentCaster;
typedef itk::IntensityWindowingImageFilter <ChannelImageType, InputComponentType>            IntensityWindowingChannelFilterType;
typedef itk::StatisticsImageFilter <ChannelImageType>                                        StatisticsChannelFilterType;
typedef itk::SubtractImageFilter< ChannelImageType, ChannelImageType, ChannelImageType >     AddConstantFilter;
typedef itk::MultiplyImageFilter< ChannelImageType, ChannelImageType, ChannelImageType >     MultiplyFilter;

// Probability Image (priors) image definitions
typedef float                                                                           ProbabilityPixelType;
typedef itk::Image< ProbabilityPixelType, 3 >                                           ProbabilityImageType;
typedef ProbabilityImageType::Pointer                                                   ProbabilityImagePointer;
typedef ProbabilityImageType::ConstPointer                                              ProbabilityImageConstPointer;
typedef itk::ImageFileReader< ProbabilityImageType >                                    ProbabilityImageReader;
typedef itk::ImageFileWriter< ProbabilityImageType >                                    ProbabilityImageWriter;

typedef itk::MaskImageFilter< ChannelImageType, ProbabilityImageType >                       MaskChannelFilterType;
typedef itk::BinaryThresholdImageFilter< InputComponentType, ProbabilityImageType >          ComponentThresholdFilter;
typedef itk::MultiplyImageFilter
		< ProbabilityImageType, ProbabilityImageType, ProbabilityImageType >             MultiplyProbabilityImage;

typedef itk::Image< unsigned char,3>                                                    ClassifiedImageType;
typedef itk::ImageFileWriter< ClassifiedImageType >                                     ClassifiedImageWriter;
typedef itk::CastImageFilter< ClassifiedImageType, ProbabilityImageType >               ClassifiedImageCaster;
typedef itk::BinaryThresholdImageFilter
		          < ProbabilityImageType, itk::Image<unsigned char, 3> >                ClassifiedImageClassFilter;


typedef itk::MultispectralFilter< InputComponentType, ClassifiedImageType >             MultispectralFilter;
typedef std::vector< InputComponentConstPointer >                                       InputVector;
typedef std::vector< ProbabilityImageConstPointer >                                     PriorsVector;

typedef itk::MultichannelParameterEstimator< InputComponentType >                       ParameterEstimator;
typedef itk::MultispectralEMFilter< InputComponentType >                                EMFilter;
typedef itk::MultispectralMRFFilter< InputComponentType >                               MRFFilter;
typedef itk::MultispectralKMeansFilter< InputComponentType, ProbabilityImageType >      KMeansFilter;
typedef itk::MultispectralMLClassifierFilter< InputComponentType,ProbabilityImageType > MLClassifierFilter;
typedef itk::IntensityWindowingImageFilter <ProbabilityImageType, ProbabilityImageType> IntensityWindowingImageFilterType;
typedef itk::ThresholdImageFilter <ProbabilityImageType>                                ThresholdImageFilterType;
typedef itk::StatisticsImageFilter <ProbabilityImageType>                               StatisticsImageFilterType;

typedef itk::Array< double >                                                            ParametersType;
typedef std::vector< ParametersType >                                                   ParametersVectorType;

enum INIT_MODES {
	NONE,
	MANUAL,
	KMEANS,
	EM,
	KMEANS_EM,
	MANUAL_KMEANS,
	MANUAL_EM
};

ParametersVectorType ReadParametersFile(std::string filename);

//typedef void ( *WriteImageFn )( const ProbabilityImageType*, std::string path );

#endif /* BRAIN_SEG_H_ */
