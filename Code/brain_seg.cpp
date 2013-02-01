// --------------------------------------------------------------------------------------
// File:          brain_seg.cpp
// Date:          Jul 14, 2011
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


#include "brain_seg.h"
#include "MBISConfig.h"

int main(int argc, char **argv) {
	ParametersVectorType initialParameters;
	std::string outPrefix, bmImageName, initFile, initMeansFile, mrfModelFile;
	std::string s_minimization = DEFAULT_MRF_ALG;
	float lambda1, lambda2;
	ProbabilityPixelType atlas_th, mask_th;
	unsigned int nClasses, emIterations, mrfIterations;
	std::vector<std::string> channels, priors, geometricTerm;
	std::string s_mode = "km+em";
	std::string outExt = "nii.gz";

	std::vector<double> dataOffsets;
	std::vector<double> dataFactors;

	bool useFile = false;
	bool skipMRF = false;
	bool usePVE = false;
	bool doOutputStats,doSegmentsOutput, doOutputMRFMaps, doOutputMRFIterations,
	     doOutputSteps, useExplicitPVE, skipNormalization, skipBias, doBiasOutput, doCorrectedOutput,
	     channelsAreMasked = false;
	unsigned int userBucketSize = 0;

	boost::unordered_map<std::string, INIT_MODES> modes_map;

	modes_map["none"] = NONE;
	modes_map["km"] = KMEANS;
	modes_map["em"] = EM;
	modes_map["km+em"] = KMEANS_EM;

	boost::unordered_map<std::string, typename MRFFilter::GCOptimizerType::MINIMIZATION_MODES > minimiz_map;
	minimiz_map["expansion"] = MRFFilter::GCOptimizerType::EXPANSION;
	minimiz_map["swap"] = MRFFilter::GCOptimizerType::SWAP;


	bpo::options_description general("General options");
	general.add_options()
		("help", "show help message")
		("out,o",bpo::value < std::string > (&outPrefix), "prefix for output files")
		("out-ext", bpo::value< std::string >( &outExt ), "output format, if supported by ITK")
		("brain-mask,x",bpo::value < std::string > (&bmImageName),"brain extracted mask")
		("channels-masked,M", bpo::bool_switch(&channelsAreMasked), "set this flag if channels are already masked")
		("channels,C", bpo::value< std::vector< std::string > >(&channels)->multitoken()->required(),"list of channels for multivariate segmentation")
		("class,n",	bpo::value<unsigned int>(&nClasses)->default_value(DEFAULT_NCLASSES), "number of tissue-type classes")
		("init,I", bpo::value< std::string >(&s_mode), "operation mode for initialization [none,km,em,km+em]")
		("init-file,f",bpo::value < std::string > (&initFile),"use file for initialization of tissue-type means & variances")
		("init-means",bpo::value < std::string > (&initMeansFile),"use file for initialization of tissue-type means")
		("segments,g", bpo::bool_switch(&doSegmentsOutput),"Output a posteriori probability for each class")
		//("pv", bpo::bool_switch(&useExplicitPVE), "Use explicit PVE Gaussian Model between tissues")
		("output-stats", bpo::bool_switch(&doOutputStats),"output classes statistics to CSV file")
		("output-steps", bpo::bool_switch(&doOutputSteps), "Output intermediate segmentation steps")
		("normalize-off,N", bpo::bool_switch(&skipNormalization), "Do not normalize input's intensity range")
		("brain-mask-th,t",bpo::value < ProbabilityPixelType > (&mask_th)->default_value(0.0), "mask probability threshold")
		("version,v", "show tool version");

	bpo::options_description kmeans_desc("Kmeans options");
	kmeans_desc.add_options()
		("bucket-size", bpo::value< unsigned int >(&userBucketSize), "bucket size for K-means operation");

	bpo::options_description em_desc("Expectation-Maximization options");
	em_desc.add_options()
		("priors,P", bpo::value< std::vector< std::string > >(&priors)->multitoken(),"add prior to list of atlases for initializing segmentation")
		("atlas-threshold",bpo::value < ProbabilityPixelType > (&atlas_th)->default_value(DEFAULT_ATLAS_TH), "threshold for atlases")
		("em-iterations", bpo::value<unsigned int>(&emIterations)->default_value(DEFAULT_MAX_ITER), "maximum number of iterations");

	bpo::options_description bias_desc("Bias estimator options");
	bias_desc.add_options()
		("bias-skip", bpo::bool_switch(&skipBias), "Do not estimate Bias field")
		("bias-output", bpo::bool_switch(&doBiasOutput), "Output estimated Bias field")
		("bias-corrected-output", bpo::bool_switch(&doCorrectedOutput), "Output corrected input images with estimated Bias field");

	bpo::options_description mrf_desc( "MRF Segmentation options" );
	mrf_desc.add_options()
		("mrf-lambda,l", bpo::value<float>(&lambda1)->default_value(DEFAULT_LAMBDA),"Regularization weighting. The higher this value is, the higher smoothing. You may use a value around 0.3-1.0 for typical brain images")
		("mrf-minimization,m", bpo::value < std::string > (&s_minimization)->default_value(DEFAULT_MRF_ALG), "Minimization algorithm")
		("mrf-energy-model,e",bpo::value < std::string > (&mrfModelFile), "Energy Model matrix, basic Pott's Model used if missing" )
		("mrf-iterations", bpo::value<unsigned int>(&mrfIterations)->default_value(DEFAULT_MRF_ITERATIONS),"MRF re-estimation iterations")
		("mrf-skip,S", bpo::bool_switch(&skipMRF),"Skip MRF step")
		("mrf-pve", bpo::bool_switch(&usePVE), "Compute PVE classes")
		("mrf-external-lambda", bpo::value<float>(&lambda2)->default_value(1.0),"External field energy weighting. The higher this value is, the higher impact from the external field")
		("mrf-external-field,E", bpo::value< std::vector< std::string > >(&geometricTerm)->multitoken(),"External field maps that manipulate the dataterm in MRF energy function")
		("mrf-output-maps", bpo::bool_switch(&doOutputMRFMaps), "Output GC probabilities and energy maps")
		("mrf-output-it", bpo::bool_switch(&doOutputMRFIterations), "Output intermediate segmentation steps");


	bpo::options_description all("Usage");
	all.add(general).add(kmeans_desc).add(bias_desc).add(em_desc).add(mrf_desc);

	bpo::variables_map vmap;
	bpo::store(bpo::command_line_parser(argc, argv).options(all).run(), vmap);


	if ( vmap.count("version")) {
			std::cout << "MBIS Segmentation tool. Version " << MBIS_VERSION_MAJOR << "." << MBIS_VERSION_MINOR << "-" << MBIS_RELEASE << "." << std::endl;
			std::cout << "--------------------------------------------------------"<< std::endl;
			std::cout << "Copyright (c) 2012, code@oscaresteban.es (Oscar Esteban)"<< std::endl; 
			std::cout << "with Signal Processing Lab 5, EPFL (LTS5-EPFL)"<< std::endl;
			std::cout << "and Biomedical Image Technology, UPM (BIT-UPM)"<< std::endl;
			std::cout << "All rights reserved."<< std::endl;

			return 0;
	}

	if ( vmap.empty() || vmap.count("help")) {
		std::cout << all << std::endl;
		return 0;
	}

	try {
		bpo::notify(vmap);
	} catch (boost::exception_detail::clone_impl<
			boost::exception_detail::error_info_injector<
					boost::program_options::required_option> > &err) {
		std::cout << "Error parsing options:" << err.what()
				<< std::endl;
		std::cout << std::endl << all << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << std::endl << "Set-up step --------------------------------------------" << std::endl;

	ProbabilityImagePointer bm;
	InputVector input;
	PriorsVector atlas;
	ClassifiedImageType::Pointer solution;
	bool useFileMask = vmap.count("brain-mask") && bfs::exists(bmImageName);
	bool useImplicitMasks = channelsAreMasked && !useFileMask;

	if( useFileMask ) {
		ProbabilityImageReader::Pointer r = ProbabilityImageReader::New();
		r->SetFileName(bmImageName);
		r->Update();

		try {
			bm = r->GetOutput();
		} catch (...) {
			std::cout << "Error reading brain mask" << std::endl;
			return EXIT_FAILURE;
		}

		StatisticsImageFilterType::Pointer calc = StatisticsImageFilterType::New();
		calc->SetInput(bm);
		calc->Update();
		ProbabilityPixelType max = calc->GetMaximum();

		if (max > 1.0 ){
			ProbabilityPixelType min = calc->GetMinimum();

			IntensityWindowingImageFilterType::Pointer intensityFilter = IntensityWindowingImageFilterType::New();
			intensityFilter->SetInput( bm );
			intensityFilter->SetWindowMaximum( max );
			intensityFilter->SetWindowMinimum( min );
			intensityFilter->SetOutputMaximum( 1.0 );
			intensityFilter->SetOutputMinimum( 0.0 );
			intensityFilter->Update();
			bm = intensityFilter->GetOutput();
		}

		std::cout << "\t* Mask: read from file " << bmImageName << ".";

		if ( mask_th != 0.0 ) {
			ThresholdImageFilterType::Pointer th = ThresholdImageFilterType::New();
			th->SetInput( bm );
			th->ThresholdBelow( mask_th );
			th->Update();
			bm = th->GetOutput();
			std::cout << " Mask Threshold = " << mask_th << ".";
		}

		std::cout << std::endl;

	} else {
		if ( useImplicitMasks ) {
			std::cout << "\t* Mask: channels are masked." << std::endl;
		} else {
			std::cout << "\t* Mask: not requested." << std::endl;
		}
	}

	std::cout << "\t* Inputs normalization is " << ((!skipNormalization)?"ON":"OFF") << std::endl;

	for (vector<string>::iterator it = channels.begin(); it != channels.end(); it++) {
		ChannelReader::Pointer r = ChannelReader::New();
		r->SetFileName( *it );
		ChannelPointer p = r->GetOutput();
		r->Update();

		ChannelPixelType max = itk::NumericTraits< ChannelPixelType >::max();
		ChannelPixelType min = 0.0;
		ChannelPixelType absMin = 0.0;

		if ( bm.IsNotNull() ) {
			ProbabilityPixelType* mBuff = bm->GetBufferPointer();
			ChannelPixelType* cBuff = r->GetOutput()->GetBufferPointer();
			size_t nPix = bm->GetLargestPossibleRegion().GetNumberOfPixels();
			std::vector< ChannelPixelType > sample;

			for( size_t i = 0; i<nPix; i++ ) {
				if ( *(mBuff+i) > 0 ) {
					sample.push_back( *(cBuff+i) );
				}
			}
			std::sort( sample.begin(), sample.end() );
			max = sample[ (size_t) ((sample.size()-1)*0.98) ];
			min = sample[ (size_t) ((sample.size()-1)*0.02) ];
			absMin = sample[0];
		} else {
			StatisticsChannelFilterType::Pointer calc = StatisticsChannelFilterType::New();
			calc->SetInput(p);
			calc->Update();
			max = calc->GetMaximum();
			min = calc->GetMinimum();
			absMin = min;
		}

		if( !skipNormalization ) {
			double factor = NORM_MAX_INTENSITY / (max - min);
			double constant = - factor * absMin;

			if ( factor!= 1 ) {
				typename MultiplyFilter::Pointer multiplier = MultiplyFilter::New();
				multiplier->SetInput( p );
				multiplier->SetConstant( factor );
				multiplier->Update();
				p = multiplier->GetOutput();
			}

			if ( constant!= 0 ) {
				typename AddConstantFilter::Pointer adder = AddConstantFilter::New();
				adder->SetInput( p );
				adder->SetConstant( - constant );
				adder->Update();
				p = adder->GetOutput();
			}

			dataOffsets.push_back( constant );
			dataFactors.push_back( factor );

			if ( bm.IsNotNull() ) {
				ChannelImageType::DirectionType dir = bm->GetDirection();
				p->SetDirection( dir );
				p->SetOrigin( bm->GetOrigin() );

				MaskChannelFilterType::Pointer m = MaskChannelFilterType::New();
				m->SetInput( p );
				m->SetMaskImage( bm );
				m->Update();
				p = m->GetOutput();
			}


			if( doOutputSteps ) {
				itk::ImageFileWriter< ChannelImageType >::Pointer wc = itk::ImageFileWriter< ChannelImageType >::New();
				wc->SetInput( p );
				std::stringstream ss;
				ss << outPrefix << "_normin_" << input.size() << "." << outExt;
				wc->SetFileName( ss.str() );
				wc->Update();
			}
		}

		InputComponentConstPointer im;

		typename ChannelToComponentCaster::Pointer c = ChannelToComponentCaster::New();
		c->SetInput(p);
		c->Update();
		im = c->GetOutput();

		input.push_back( im );

		std::cout << "\t* Channel [" << std::setw(2) << input.size() << "] read: " << (*it) << std::endl;
	}

	if ( useImplicitMasks ) {
		bm = ProbabilityImageType::New();
		bm->SetRegions( input[0]->GetLargestPossibleRegion() );
		bm->CopyInformation( input[0] );
		bm->Allocate();
		bm->FillBuffer( 1.0 );
		bm->Update();

		ProbabilityPixelType* buffer = bm->GetBufferPointer();
		size_t numberOfPixels = bm->GetLargestPossibleRegion().GetNumberOfPixels();

		for ( size_t i = 0; i < input.size(); i++) {
			const PixelType* buffer2 = input[i]->GetBufferPointer();

			for( size_t offset = 0; offset < numberOfPixels; offset++ ) {
				if( *( buffer + offset ) > 0.0 ) {
					*( buffer + offset ) = PixelType ( *( buffer2 + offset ) > 1e-5);
				}
			}
		}


		if( doOutputSteps ) {
			itk::ImageFileWriter< ProbabilityImageType >::Pointer wc = itk::ImageFileWriter< ProbabilityImageType >::New();
			wc->SetInput( bm );
			std::stringstream ss;
			ss << outPrefix << "_mask." << outExt;
			wc->SetFileName( ss.str() );
			wc->Update();
		}
	}


	unsigned int init_mode = modes_map[s_mode];
	unsigned int nComponents = input.size();
	unsigned int nElements = nComponents * (1+nComponents);

	if( priors.size() > 0 ) {
		if (init_mode== KMEANS ) init_mode = MANUAL;
		if (init_mode == KMEANS_EM ) init_mode = EM;

		for (vector<string>::iterator it = priors.begin(); it != priors.end(); it++) {
			ProbabilityImageReader::Pointer r = ProbabilityImageReader::New();
			r->SetFileName( *it );
			ProbabilityImageType::ConstPointer p = r->GetOutput();
			r->Update();

			StatisticsImageFilterType::Pointer calc = StatisticsImageFilterType::New();
			calc->SetInput(p);
			calc->Update();
			ProbabilityPixelType max = calc->GetMaximum();

			if (max > 1.0 ){
				ProbabilityPixelType min = calc->GetMinimum();

				IntensityWindowingImageFilterType::Pointer intensityFilter = IntensityWindowingImageFilterType::New();
				intensityFilter->SetInput( p );
				intensityFilter->SetWindowMaximum( max );
				intensityFilter->SetWindowMinimum( min );
				intensityFilter->SetOutputMaximum( 1.0 );
				intensityFilter->SetOutputMinimum( 0.0 );
				intensityFilter->Update();
				p = intensityFilter->GetOutput();
			}
			atlas.push_back(p);

			std::cout << "\t* Prior [" << std::setw(2) << atlas.size() << "] read: " << (*it) << std::endl;

			ParameterEstimator::Pointer est = ParameterEstimator::New();
			est->SetInputVector( input );
			est->SetPrior( p );
			if ( bm.IsNotNull() ) {
				est->SetMaskImage( bm );
			}

			est->Update();
			initialParameters.push_back( est->GetOutputParameters() );
		}
	}

	if ( bfs::exists(initFile)) {
		if (init_mode== KMEANS ) init_mode = MANUAL;
		if (init_mode == KMEANS_EM ) init_mode = EM;

		std::cout << "\t* Parsing tissue parameters file: " << initFile << std::endl;
		initialParameters = ReadParametersFile(initFile);

		for ( ParametersVectorType::iterator it = initialParameters.begin(); it!=initialParameters.end(); it++) {
			size_t n = (*it).size();
			if ( n != nElements ) {
				std::cerr << "Parameters file is incorrect or badly interpreted" << std::endl;
			}

			if ( !skipNormalization ) {
				for( size_t i = 0; i<nComponents; i++ ) {
					(*it)[i] *= dataFactors[i];
					(*it)[i] += dataOffsets[i];
				}

				for( size_t i = nComponents; i<n; i++ ) {
					(*it)[i] = dataFactors[i%nComponents] * (*it)[i];
				}
			}
		}

		useFile = true;

		std::cout << "\t* Manual Parameters: " << std::endl;
		for (size_t i = 0; i<initialParameters.size(); i++)
			std::cout << "\t\t" << initialParameters[i] << std::endl;
	}

	if ( bfs::exists(initMeansFile)) {
		std::cout << "\t* Parsing tissue means file: " << initMeansFile << std::endl;

		initialParameters = ReadParametersFile( initMeansFile );

		for ( ParametersVectorType::iterator it = initialParameters.begin(); it!=initialParameters.end(); it++) {
			size_t n = (*it).size();

			if ( n != nComponents ) {
				std::cerr << "Means file is incorrect or badly interpreted" << std::endl;
			}


			if ( !skipNormalization ) {
				for( size_t i = 0; i<n; i++ ) {
					(*it)[i] *= dataFactors[i];
					(*it)[i] += dataOffsets[i];
				}
			}
		}

		useFile = true;

		std::cout << "\t* Manual Parameters: " << std::endl;
		for (size_t i = 0; i<initialParameters.size(); i++)
			std::cout << "\t\t" << initialParameters[i] << std::endl;
	}

	EMFilter::Pointer em_filter;
	KMeansFilter::Pointer kmeans;


	if (init_mode == KMEANS || init_mode == KMEANS_EM ) {
		std::cout << std::endl << "Kmeans step --------------------------------------------" << std::endl;
		kmeans = KMeansFilter::New();
		kmeans->SetInputVector( input );
		if ( bm.IsNotNull() ) kmeans->SetMaskImage( bm );
		kmeans->SetNumberOfClasses(nClasses);
		if(vmap.count("bucket-size")) {
			kmeans->SetUserBucketSize(userBucketSize);
		}

		if( useFile ) {
			kmeans->SetInitialParameters( initialParameters );
		}

		kmeans->Compute();
		initialParameters = kmeans->GetOutputParameters();

		if (doOutputStats) {
			stringstream s;
			s << outPrefix << "_stats_kmeans.csv";
			std::ofstream outputCSVfile ( s.str().c_str() );
			for ( unsigned int i = 0; i < initialParameters.size(); i++)
				outputCSVfile << initialParameters[i] << "\n";
			outputCSVfile.close();
		}
		kmeans->Update();
	}

	// Check for sanity initial parameters
	if( initialParameters.size()!=nClasses ) {
		std::cerr << "Error: initial parameters size (" << initialParameters.size() << ") doesn't match number of classes (" << (int) nClasses << std::endl;
	} else if ( initialParameters[0].size() != nElements ) {
		typename MLClassifierFilter::Pointer ml_classifier = MLClassifierFilter::New();
		ml_classifier->SetInputVector( input );
		if ( bm.IsNotNull() ) ml_classifier->SetMaskImage(bm);
		ml_classifier->SetParameters( initialParameters );
		ml_classifier->Update();
		solution = ml_classifier->GetOutput();

		if ( doOutputSteps ) {
			ClassifiedImageWriter::Pointer w = ClassifiedImageWriter::New();
			w->SetInput( solution );
			w->SetFileName(outPrefix + "_means." + outExt );
			w->Update();
		}

		ProbabilityImagePointer unoImage = ProbabilityImageType::New();
		unoImage->SetRegions( input[0]->GetLargestPossibleRegion() );
		unoImage->CopyInformation( input[0] );
		unoImage->Allocate();
		unoImage->FillBuffer( 1.0 );
		unoImage->Update();

		// Initialize Covariance matrix --------------------------------------
		ParameterEstimator::Pointer est = ParameterEstimator::New();
		est->SetInputVector( input );

		typedef itk::LabelImageToLabelMapFilter< ClassifiedImageType > ClassifiedToLabelMapFilter;
		typedef typename ClassifiedToLabelMapFilter::OutputImageType LabelMap;
		typedef itk::LabelMapMaskImageFilter< LabelMap, ProbabilityImageType > LabelMapMaskFilter;

		unsigned int maskOffset = (unsigned int) bm.IsNotNull();

	    for ( typename ClassifiedImageType::PixelType i = 0 ; i < nClasses ; i++ ) {
			typename ClassifiedToLabelMapFilter::Pointer lfilter = ClassifiedToLabelMapFilter::New();
			lfilter->SetInput( solution );
			lfilter->Update();

			typename LabelMapMaskFilter::Pointer lmask = LabelMapMaskFilter::New();
			lmask->SetInput1( lfilter->GetOutput() );
			lmask->SetInput2( unoImage );
			lmask->SetLabel( i + maskOffset );
			lmask->Update();

			est->SetPrior( lmask->GetOutput());
			est->Update();
			initialParameters[i] = est->GetOutputParameters();
	    }
		// End initialize covariance matrix -------------------------------------------
	}

	if (init_mode == EM  || init_mode == KMEANS_EM ) {
		std::cout << std::endl << "E-M Step -----------------------------------------------" << std::endl;
		em_filter = EMFilter::New();
		em_filter->SetMaskImage( bm );
		em_filter->SetNumberOfClasses( nClasses );
		em_filter->SetMaximumIteration( emIterations );
		em_filter->SetMaxBiasEstimationIterations( 5 );
		em_filter->SetInputVector( input );
		em_filter->SetInitialParameters( initialParameters );
		em_filter->SetUseExplicitPVModel( useExplicitPVE );
		em_filter->SetUseBiasCorrection( !skipBias );

		if ( atlas.size() != 0 ) {
			em_filter->SetPriors( atlas );
			em_filter->SetUsePriorProbabilityImages( true );
		}

		em_filter->Update();

		// TODO change to GetModelParameters()
		initialParameters = em_filter->GetInitialParameters();
		solution = em_filter->GetOutput();

		if ( doOutputSteps || skipMRF ) {
			ClassifiedImageWriter::Pointer w = ClassifiedImageWriter::New();
			w->SetInput( solution );
			w->SetFileName(outPrefix + "_em." + outExt );
			w->Update();
		}

		//
		// Output file with initial stats if switch is true
		//
		if (doOutputStats) {
			stringstream s;
			s << outPrefix << "_stats_em.csv";
			std::ofstream outputCSVfile ( s.str().c_str() );
			for ( unsigned int i = 0; i < initialParameters.size(); i++)
				outputCSVfile << initialParameters[i] << "\n";
			outputCSVfile.close();
		}

		if( em_filter->GetUseBiasCorrection() ) {
			typedef typename EMFilter::EstimatorType::InputVectorImageType  VectorImage;
			typedef typename itk::ImageFileWriter< ChannelImageType > ChannelWriter;
			typedef itk::VectorIndexSelectionCastImageFilter< VectorImage, ChannelImageType> SelectorType;
			typename VectorImage::ConstPointer vec = em_filter->GetEstimator()->GetCorrectedInput();
			for( int channel = 0; channel< input.size(); channel++ ) {
				typename SelectorType::Pointer channelSelector = SelectorType::New();
				channelSelector->SetInput( vec );
				channelSelector->SetIndex( channel );
				channelSelector->Update();
				channelSelector->GetOutput()->SetRegions( vec->GetRequestedRegion() );
				input[channel] = channelSelector->GetOutput();

				if ( doCorrectedOutput ) {
					char name[50];
					sprintf(name, "_corrected_ch%02d.%s" , channel, outExt.c_str() );
					typename ChannelWriter::Pointer chW = ChannelWriter::New();
					chW->SetInput(input[channel]);
					chW->SetFileName( outPrefix + name );
					chW->Update();
				}
			}

			if( doBiasOutput ) {
				typedef typename EMFilter::EstimatorType::InputVectorImageType  VectorImage;
				typedef typename itk::ImageFileWriter< ChannelImageType > ChannelWriter;
				typename VectorImage::ConstPointer bias = em_filter->GetEstimator()->GetCurrentBias();
				for( int channel = 0; channel< input.size(); channel++ ) {
					typename SelectorType::Pointer channelSelector = SelectorType::New();
					channelSelector->SetInput( bias );
					channelSelector->SetIndex( channel );
					channelSelector->Update();
					channelSelector->GetOutput()->SetRegions( bias->GetRequestedRegion() );

					char name[50];
					sprintf(name, "_bias_ch%02d.%s" , channel, outExt.c_str() );
					typename ChannelWriter::Pointer chW = ChannelWriter::New();
					chW->SetInput(channelSelector->GetOutput());
					chW->SetFileName( outPrefix + name );
					chW->Update();
				}
			}

		}
	}

	//
	// MRF Segmentation
	//

	if( !skipMRF ) {
		std::cout << std::endl << "MRF Step -----------------------------------------------" << std::endl;
		MRFFilter::Pointer mrf = MRFFilter::New();
		mrf->SetUseExplicitPVModel( useExplicitPVE );
		mrf->SetNumberOfClasses( nClasses );
		mrf->SetMaskImage( bm );
		mrf->SetInputVector( input );
		mrf->SetLambda( lambda1 );
		mrf->SetExternalLambda( lambda2 );
		mrf->SetUseOutputPriors( doSegmentsOutput );
		mrf->SetInitialParameters( initialParameters );
		mrf->SetMRFIterations( mrfIterations );
		mrf->SetMinimizationMode( minimiz_map[s_minimization] );
		mrf->SetOutputPrefix( outPrefix );

		if( doOutputMRFIterations ) {
			typedef itk::SimpleMemberCommand< MRFFilter >  ObserverType;
			ObserverType::Pointer obs = ObserverType::New();
			obs->SetCallbackFunction( mrf, &MRFFilter::WriteIteration );
			mrf->AddObserver( itk::IterationEvent(), obs );
		}


		if( doOutputMRFMaps ) {
			typedef itk::SimpleMemberCommand< MRFFilter >  ObserverType;
			ObserverType::Pointer obs = ObserverType::New();
			obs->SetCallbackFunction( mrf, &MRFFilter::WriteMaps );
			mrf->AddObserver( itk::IterationEvent(), obs );
		}

		if ( doOutputSteps && init_mode == NONE ) {
			EMFilter::Pointer em_classifier = EMFilter::New();
			em_classifier->SetMaskImage( bm );
			em_classifier->SetNumberOfClasses( nClasses );
			em_classifier->SetInputVector( input );
			em_classifier->SetInitialParameters( initialParameters );
			em_classifier->SetUseExplicitPVModel( useExplicitPVE );
			em_classifier->SetUseOnlyClassify(true);
			em_classifier->Update();
			ClassifiedImageWriter::Pointer w = ClassifiedImageWriter::New();
			w->SetInput( em_classifier->GetOutput() );
			w->SetFileName(outPrefix + "_mrf_initialization." + outExt );
			w->Update();

		}

		if( bfs::exists( mrfModelFile ) ) {
			std::cout << "\t* Loading Energy Model Matrix: " << mrfModelFile << std::endl;
			mrf->SetMatrixFile( mrfModelFile );
		}

		for (vector<string>::iterator it = geometricTerm.begin(); it != geometricTerm.end(); it++) {
			if( bfs::exists(*it) ) {
				size_t label = (it - geometricTerm.begin()) + 1;

				ProbabilityImageReader::Pointer r = ProbabilityImageReader::New();
				r->SetFileName( *it );
				ProbabilityImageType::Pointer p = r->GetOutput();
				r->Update();

				mrf->SetSpatialEnergyMap(label, p );


				std::cout << "\t* Geometrical constraint [" << std::setw(2) << label << "] read: " << (*it) << std::endl;
			}
		}

		mrf->Update();
		solution = mrf->GetOutput();

		ClassifiedImageWriter::Pointer w2 = ClassifiedImageWriter::New();
		w2->SetInput( solution );
		w2->SetFileName(outPrefix + "_mrf." + outExt );
		w2->Update();

		if ( doSegmentsOutput ) {
			for ( unsigned int i = 1; i<=nClasses; i++ ) {
				ProbabilityImageWriter::Pointer wProb = ProbabilityImageWriter::New();
				wProb->SetInput( mrf->GetPosteriorProbabilityImage(i) );
				char name[50];
				sprintf(name, "_mrf_seg%02d.%s" , i, outExt.c_str() );
				wProb->SetFileName( outPrefix + name );
				wProb->Update();
			}

		}

		if (doOutputStats) {
			initialParameters = mrf->GetInitialParameters();
			stringstream s;
			s << outPrefix << "_stats_final.csv";
			std::ofstream outputCSVfile ( s.str().c_str() );
			for ( unsigned int i = 0; i < initialParameters.size(); i++)
				outputCSVfile << initialParameters[i] << "\n";
			outputCSVfile.close();
		}
	}


	return EXIT_SUCCESS;
}

ParametersVectorType ReadParametersFile(std::string filename) {
	ParametersVectorType result;

	std::ifstream file(filename.c_str());
	std::string line;

	typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
	boost::char_separator<char> sep(",;| []");

	unsigned int classId = 0;

	while ( getline(file, line)) {
		size_t pos = line.find("#");

		if( pos!=0 ) {
			line = line.substr(0,pos);

			tokenizer toker(line, sep);
			tokenizer::iterator iterator = toker.begin();
			std::vector<float> classStats;

			while (iterator != toker.end()) {
				float val = boost::lexical_cast<float>(*iterator++);
				classStats.push_back(val);
			}

			ParametersType statsArray(classStats.size());
			std::copy(classStats.begin(), classStats.end(), statsArray.begin());
			result.push_back(statsArray);
			classId++;
		}
	}

	file.close();

	return result;
}
