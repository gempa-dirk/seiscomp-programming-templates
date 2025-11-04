/***************************************************************************
 * Copyright (C) Jan Becker, gempa GmbH                                    *
 * All rights reserved.                                                    *
 * Contact: jabe@gempa.de                                                  *
 *                                                                         *
 * GNU Affero General Public License Usage                                 *
 * This file may be used under the terms of the GNU Affero                 *
 * Public License version 3.0 as published by the Free Software Foundation *
 * and appearing in the file LICENSE included in the packaging of this     *
 * file. Please review the following information to ensure the GNU Affero  *
 * Public License version 3.0 requirements will be met:                    *
 * https://www.gnu.org/licenses/agpl-3.0.html.                             *
 ***************************************************************************/


// Defines the logging component which can be shown when using
// --print-component 1 with each SeisComP application.
#define SEISCOMP_COMPONENT PGA

// Include the respective headers
// - Logging: to use the SEISCOMP_* logging macros
#include <seiscomp/logging/log.h>
// - Plugin: to enable building a plugin
#include <seiscomp/core/plugin.h>
// - Filtering
#include <seiscomp/math/filter/chainfilter.h>
#include <seiscomp/math/filter/iirdifferentiate.h>
#include <seiscomp/math/filter/iirintegrate.h>
#include <seiscomp/math/filter/rmhp.h>
// - AmplitudeProcessor interface required to derive from it
#include <seiscomp/processing/amplitudeprocessor.h>
// - For channel data combiners
#include <seiscomp/processing/operator/l2norm.h>
#include <seiscomp/processing/operator/ncomps.h>


// Everything is implemented in a private namespace to not export any symbols to keep
// the global symbol space clean when loading the plugin.
namespace {


// Use some standard namespaces to write more compact code.
using namespace std;
using namespace Seiscomp;
using namespace Seiscomp::Processing;


#define AMPLITUDE_TYPE "template_pga"


class PGAProcessor : public AmplitudeProcessor {
	// ----------------------------------------------------------------------
	//  X'truction
	// ----------------------------------------------------------------------
	public:
		/**
		 * @brief The constructor.
		 * The amplitude type can be an arbitrary name but should not conflict
		 * with any other existing name. It should be the same as the registered
		 * interface name, see the REGISTER_AMPLITUDEPROCESSOR statement at the
		 * end of this file.
		 */
		PGAProcessor()
		: AmplitudeProcessor(AMPLITUDE_TYPE) {
			// Configure the relative time window of the data with respect to
			// the set trigger time.
			setNoiseStart(-10);
			setNoiseEnd(-2);
			setSignalStart(-2);
			setSignalEnd("max(150, R / 3.5)");

			// Tell the amplitude processor to feed data for both horizontal
			// components.
			setDataComponents(Horizontal);
			setTargetComponent(FirstHorizontalComponent);
		}


	// ----------------------------------------------------------------------
	//  Public AmplitudeProcessor interface
	// ----------------------------------------------------------------------
	public:
		bool setup(const Settings &settings) override {
			// Reset operator and filter
			setOperator(nullptr);
			setFilter(nullptr);
			clear(State::FilterCreated);

			// Call base class implementation
			if ( !AmplitudeProcessor::setup(settings) ) {
				// If the setup of the base class fails, then it does not make sense
				// to continue.
				return false;
			}

			// Check the horizontal components for valid gains
			for ( int i = FirstHorizontal; i <= SecondHorizontal; ++i ) {
				if ( _streamConfig[i].code().empty() ) {
					SEISCOMP_ERROR("Component[%d] code is empty", i);
					setStatus(Error, i);
					return false;
				}

				if ( _streamConfig[i].gain == 0.0 ) {
					SEISCOMP_ERROR("Component[%d] gain is missing (actually zero)", i);
					setStatus(MissingGain, i);
					return false;
				}
			}

			if ( _streamConfig[FirstHorizontal].gainUnit != _streamConfig[SecondHorizontal].gainUnit ) {
				SEISCOMP_ERROR("Both components do not have the same gain unit: %s != %s",
				               _streamConfig[FirstHorizontal].gainUnit,
				               _streamConfig[SecondHorizontal].gainUnit);
				setStatus(ConfigurationError, 1);
				return false;
			}

			string preFilter, postFilter;

			try { preFilter = settings.getString("amplitudes." + type() + ".preFilter"); }
			catch ( ... ) {}

			try { postFilter = settings.getString("amplitudes." + type() + ".filter"); }
			catch ( ... ) {}

			SEISCOMP_DEBUG("  + pre-filter = %s", preFilter);
			SEISCOMP_DEBUG("  + filter = %s", postFilter);

			using OpWrapper = Operator::StreamConfigWrapper<double, 2, Operator::L2Norm>;

			if ( !preFilter.empty() ) {
				// Create a filter instance from the provided string.
				string error;
				auto filter = Filter::Create(preFilter, &error);
				if ( !filter ) {
					// If the string is wrong
					SEISCOMP_ERROR("Failed to create pre-filter: %s: %s", preFilter, error);
					setStatus(ConfigurationError, 2);
					return false;
				}

				// Create a waveform operator that combines the two horizontal channels
				// and computes L2 of each horizontal component filtered sample.
				using FilterL2Norm = Operator::FilterWrapper<double, 2, OpWrapper>;
				setOperator(
					new NCompsOperator<double, 2, FilterL2Norm>(
						FilterL2Norm(
							filter, OpWrapper(
								_streamConfig + FirstHorizontal,
								Operator::L2Norm<double, 2>()
							)
						)
					)
				);
			}
			else {
				// Create a waveform operator that combines the two horizontal channels
				// and computes L2 of each horizontal component sample.
				setOperator(
					new NCompsOperator<double, 2, OpWrapper>(
						OpWrapper(
							_streamConfig + FirstHorizontal,
							Operator::L2Norm<double, 2>()
						)
					)
				);
			}


			if ( !postFilter.empty() ) {
				string error;
				auto filter = Filter::Create(postFilter, &error);
				if ( !filter ) {
					// If the string is wrong
					SEISCOMP_ERROR("Failed to create filter: %s: %s", preFilter, error);
					setStatus(ConfigurationError, 3);
					return false;
				}

				setFilter(filter);
			}


			if ( !_enableResponses ) {
				// If responses are not enabled then we need to convert the data to
				// the required unit in time domain with IIR filters. So we have to
				// form a new filter which also takes the optional current filter into
				// account.

				SignalUnit unit;
				// Valid value already checked in setup()
				unit.fromString(_streamConfig[FirstHorizontal].gainUnit.c_str());

				Filter *filter{nullptr};
				Math::Filtering::ChainFilter<double> *chain{nullptr};

				if ( unit == Meter ) {
					SEISCOMP_DEBUG("Add double derivation of data for PGA computation");
					chain = new Math::Filtering::ChainFilter<double>;
					chain->add(new Math::Filtering::IIRDifferentiate<double>);
					chain->add(new Math::Filtering::IIRDifferentiate<double>);
				}
				else if ( unit == MeterPerSecond ) {
					SEISCOMP_DEBUG("Add single derivation of data for PGA computation");
					filter = new Math::Filtering::IIRDifferentiate<double>;
				}

				if ( chain || filter ) {
					if ( _stream.filter ) {
						if ( !chain ) {
							chain = new Math::Filtering::ChainFilter<double>;
						}

						if ( filter ) {
							chain->add(filter);
						}

						chain->add(_stream.filter);
						filter = nullptr;
					}

					_stream.filter = filter ? filter : chain;
				}

				// Set the flag that the filter chain for data conversion has been
				// created already so that the base AmplitudeProcessor does not
				// create the filter again in initFilter().
				set(State::FilterCreated);
			}

			return true;
		}


		void prepareData(DoubleArray &data) override {
			Sensor *sensor = _streamConfig[FirstHorizontal].sensor();

			// When using full responses then all information needs to be set up
			// correctly otherwise an error is set
			if ( _enableResponses && !check(State::ResponseApplied) ) {
				if ( !sensor ) {
					setStatus(MissingResponse, 1);
					return;
				}

				if ( !sensor->response() ) {
					setStatus(MissingResponse, 2);
					return;
				}

				// If the unit cannot be converted into the internal
				// enum (what basically means "unknown") then the deconvolution
				// cannot be correctly. We do not want to assume a unit here
				// to prevent computation errors in case of bad configuration.
				SignalUnit unit;
				if ( !unit.fromString(_streamConfig[FirstHorizontal].gainUnit.c_str()) ) {
					// Invalid unit string
					setStatus(IncompatibleUnit, 2);
					return;
				}

				int intSteps = 0;
				switch ( unit ) {
					case Meter:
						intSteps = -2;
						break;
					case MeterPerSecond:
						intSteps = -1;
						break;
					case MeterPerSecondSquared:
						break;
					default:
						setStatus(IncompatibleUnit, 1);
						return;
				}

				set(State::ResponseApplied);

				if ( !deconvolveData(sensor->response(), _data, intSteps) ) {
					setStatus(DeconvolutionFailed, 0);
					return;
				}
			}
		}


		bool feed(const Record *rec) override {
			if ( !getOperator() ) {
				SEISCOMP_ERROR("No operator set, has setup() been called?");
				return false;
			}

			return AmplitudeProcessor::feed(rec);
		}


		//! See Seiscomp::Processing::AmplitudeProcessor::computeAmplitude for
		//! more documentation of this function. It actually computes the
		//! amplitude.
		bool computeAmplitude(const DoubleArray &data,
		                      size_t i1, size_t i2,
		                      size_t si1, size_t si2,
		                      double offset,
		                      AmplitudeIndex *dt,
		                      AmplitudeValue *amplitude,
		                      double *period, double *snr) override {
			// Data is in acceleration: m/s**2
			*period = -1;
			*snr = -1;
			dt->index = find_absmax(data.size(), data.typedData(), si1, si2, offset);
			amplitude->value = abs(data[dt->index] - offset);

			if ( *_noiseAmplitude == 0. ) {
				*snr = -1;
			}
			else {
				*snr = amplitude->value / *_noiseAmplitude;
			}

			if ( *snr < _config.snrMin ) {
				setStatus(LowSNR, *snr);
				return false;
			}

			return true;
		}
};


}


// This defines the entry point after loading this library dynamically. That is
// mandatory if a plugin should be loaded with a SeisComP application.
ADD_SC_PLUGIN(
	"Amplitude PGA plugin template, it just computes the PGA.",
	"Jan Becker, gempa GmbH",
	0, 0, 1
)


// Bind the class PGAProcessor to the name "template_pga". This allows later to
// instantiate this class via the amplitude name. This is different to the generic
// class factory as it registeres a concrete interface, here, the AmplitudeProcessor.
REGISTER_AMPLITUDEPROCESSOR(PGAProcessor, AMPLITUDE_TYPE);
