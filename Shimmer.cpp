//-------------------------------------------------------------------------------------------------------
//  Shimmer.cpp
//  BranchReverb plugin based on the Freeverb scheme
//  Customized by Renato Volpe on 07/02/2022.
//
//-------------------------------------------------------------------------------------------------------

#include "Shimmer.h"
#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <math.h>
#include "utils.h"

/*--------------------------------------------------------------------*/
// Plugin constants
#define NUM_PRESETS 5
/*--------------------------------------------------------------------*/

/*--------------------------------------------------------------------*/
// FDN constants
#define DEFAULT_NUMBER_OF_INTERNAL_CHANNELS_FDN	16
#define DIFFUSER_DELAY_BUFFER_SIZE_MS 2000.0
#define FEEDBACK_DELAY_BUFFER_SIZE_MS 2000.0
#define MAX_REVERB_DECAY_IN_SECONDS 30.0
#define MIN_FEEDBACK_DELAY_LENGTH 100.0
#define NUMBER_OF_DIFFUSION_STEPS 5
#define EARLYREFL_DELAY_DISTRIBUTION DelayDistribution::RandomInRange
#define DIFFUSER_DELAY_DISTRIBUTION DelayDistribution::RandomInRange
#define FEEDBACK_DELAY_DISTRIBUTION DelayDistribution::RandomInRange
#define OUTPUT_LPF_TYPE LPFilterType::Shelving
#define DAMPING_LPF_TYPE LPFilterType::Vicanek
#define OUTPUT_HPF_TYPE HPFilterType::Shelving
#define DIFFUSION_LOGIC DiffuserDelayLogic::Doubled
#define MIN_DAMPING_FREQUENCY 200.0
#define MAX_DAMPING_FREQUENCY 20000.0
#define MAX_MOD_RATE 5.0
#define LPF_FILTER_MIN_FREQ 100.0
#define HPF_FILTER_MIN_FREQ 40.0
#define LPF_FILTER_MAX_FREQ 20000.0
#define HPF_FILTER_MAX_FREQ 7000.0
#define BRANCH_REVERB_DECAY 6.0
const float MIN_DAMPING_FREQUENCY_LOG = log(MIN_DAMPING_FREQUENCY);
const float MAX_DAMPING_FREQUENCY_LOG = log(MAX_DAMPING_FREQUENCY);
const float LPF_FILTER_MAX_FREQ_LOG = log(LPF_FILTER_MAX_FREQ);
const float LPF_FILTER_MIN_FREQ_LOG = log(LPF_FILTER_MIN_FREQ);
/*--------------------------------------------------------------------*/

/*--------------------------------------------------------------------*/
// Pitch shifter constants
#define NUM_OF_PITCH_INTERVALS_ALLOWED 10
const float DELTA_PARAMETER_BETWEEN_INTERVALS = 1.0/NUM_OF_PITCH_INTERVALS_ALLOWED;
char* INTERVALS_NAMES_STRING[NUM_OF_PITCH_INTERVALS_ALLOWED] = { "2nd Maj", "3rd Min", "3rd Maj", "4th Per", "5th Per", "6th Maj", "7th Maj", "1st Oct", "1 Oct+5", "1+2 Oct"};
const float INTERVALS_IN_SEMITONES_PITCH1[NUM_OF_PITCH_INTERVALS_ALLOWED] = {2.0, 3.0, 4.0, 5.0, 7.0, 9.0, 11.0, 12.0, 12.0, 12.0};
const float INTERVALS_IN_SEMITONES_PITCH2[NUM_OF_PITCH_INTERVALS_ALLOWED] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0,  0.0, 19.0, 24.0};
/*--------------------------------------------------------------------*/


/*--------------------------------------------------------------------*/
// Mix dry and wet signals with equal power distribution
void Shimmer::updateMix() {
    _wet = sin(shim_mix * M_PI * 0.5);
    _dry = cos(shim_mix * M_PI * 0.5);
}

// If two pitch shifters are giving output, sum them with half gain each
void Shimmer::updateMixPitchShifters(float pitch2) {
    _mixP1 = 1.0;
    _mixP2 = 0.0;
    if (pitch2 != 0.0) {
        _mixP1 = 0.5;
        _mixP2 = 0.5;
    }
}

/*--------------------------------------------------------------------*/
// Shimmer class constructor
Shimmer::Shimmer(audioMasterCallback audioMaster)
    : AudioEffectX(audioMaster, NUM_PRESETS, Param_Count) // n program, n parameters
{
    setNumInputs(2);		// stereo in
    setNumOutputs(2);		// stereo out
    setUniqueID('Fox');	    // identify    
    InitPlugin();
}
/*--------------------------------------------------------------------*/

/*--------------------------------------------------------------------*/
AudioEffect* createEffectInstance(audioMasterCallback audioMaster)
{
    return new Shimmer(audioMaster);
}
/*--------------------------------------------------------------------*/

/*--------------------------------------------------------------------*/
// Initialize all the objects and parameters
void Shimmer::InitPlugin()
{
    // get current sample rate
    int sampleRate = getSampleRate();

    /*.......................................*/
    // Create FDN Branch Reverb
    BranchReverb = new FDN(1, DEFAULT_NUMBER_OF_INTERNAL_CHANNELS_FDN, 2, NUMBER_OF_DIFFUSION_STEPS, 1);

    // Initialize objects (allocate delay lines)
    BranchReverb->initialize(DIFFUSER_DELAY_BUFFER_SIZE_MS, FEEDBACK_DELAY_BUFFER_SIZE_MS, sampleRate);
    
    // output mixing mode
    BranchReverb->setMixMode(MixMode::First);

    // filters' types
    BranchReverb->setDampingType(DAMPING_LPF_TYPE);
    BranchReverb->setLowPassType(OUTPUT_LPF_TYPE);
    BranchReverb->setHighPassType(OUTPUT_HPF_TYPE);

    /*.......................................*/
    // Create FDN Master Reverb
    MasterReverb = new FDN(2, DEFAULT_NUMBER_OF_INTERNAL_CHANNELS_FDN, 2, NUMBER_OF_DIFFUSION_STEPS, 1);

    // Initialize objects (allocate delay lines)
    MasterReverb->initialize(DIFFUSER_DELAY_BUFFER_SIZE_MS, FEEDBACK_DELAY_BUFFER_SIZE_MS, sampleRate);

    // output mixing mode
    MasterReverb->setMixMode(MixMode::First);

    // filters' types
    MasterReverb->setDampingType(DAMPING_LPF_TYPE);
    MasterReverb->setLowPassType(OUTPUT_LPF_TYPE);
    MasterReverb->setHighPassType(OUTPUT_HPF_TYPE);

    /*.......................................*/
    // init PSMVocoder
    PitchShift_1oct = new PSMVocoder();
    PitchShift_2oct = new PSMVocoder();

    // set sample rate
    PitchShift_1oct->reset((double)sampleRate);
    PitchShift_2oct->reset((double)sampleRate);

    // set phase locking and peak tracking
    PSMVocoderParameters params1 = PitchShift_1oct->getParameters();
    PSMVocoderParameters params2 = PitchShift_2oct->getParameters();

    params1.enablePeakPhaseLocking = true;
    params2.enablePeakPhaseLocking = true;
    params1.enablePeakTracking = true;
    params2.enablePeakTracking = true;

    PitchShift_1oct->setParameters(params1);
    PitchShift_2oct->setParameters(params2);

    /*.......................................*/
    // initialize reverb plug-in parameters
    InitPresets();  
    
    /*.......................................*/
    // set fixed parameters for BranchReverb
    BranchReverb->setLowPassFrequency(LPF_FILTER_MAX_FREQ);
    BranchReverb->setHighPassFrequency(HPF_FILTER_MIN_FREQ);

    // Modulation
    BranchReverb->setModDepth(0.0);
    BranchReverb->setModRate(0.0);

    // stereo spread
    BranchReverb->setStereoSpread(0.5);          
    /*.......................................*/    
}
/*--------------------------------------------------------------------*/

/*--------------------------------------------------------------------*/
// replace the "setSampleRate" method with user-defined one
void Shimmer::setSampleRate(float sampleRate)
{
    // Call AudioEffect "setSampleRate" method
    AudioEffect::setSampleRate(sampleRate);

    // Call setSampleRate on every needed module
    BranchReverb->setSampleRate(sampleRate);
    MasterReverb->setSampleRate(sampleRate);
    PitchShift_1oct->reset(sampleRate);
    PitchShift_2oct->reset(sampleRate);
}
/*--------------------------------------------------------------------*/

/* ------------------------------------------------------------------------------------------------------------
  -------------------------------------  PROCESS REPLACING  ---------------------------------------------------
  ------------------------------------------------------------------------------------------------------------ */
void Shimmer::processReplacing(float** inputs, float** outputs, VstInt32 sampleFrames)
{
    // Extract input and output buffers
    float* inL = inputs[0]; // buffer input left
    float* inR = inputs[1]; // buffer input right

    float* outL = outputs[0]; // buffer output left
    float* outR = outputs[1]; // buffer output right

    // write input to file
    //string pre = "test_input.txt";
    //WriteBufferToFile(inputs, sampleFrames, pre); 
    
    // Cycle over the sample frames number
    for (int i = 0; i < sampleFrames; i++) {

        // Create tmp arrays for processing
        float pitchIn = (inL[i] + inR[i]) / 2;
        double pitchOut1 = 0.0;
        double pitchOut2 = 0.0;
        float pitchOut;
        float bran_rev_out[2] = { 0.0, 0.0 };
        float mast_rev_out[2] = { 0.0, 0.0 };
        float mast_rev_in[2];

        // --- Pitch Shifting        
        // Process pitch shifting 1 octave
        pitchOut1 = PitchShift_1oct->processAudioSample(pitchIn);

        // Process pitch shifting 2 octaves
        pitchOut2 = PitchShift_2oct->processAudioSample(pitchIn);

        // Sum outputs
        pitchOut = _mixP1 * (float)pitchOut1 + _mixP2 * (float)pitchOut2;
              
        // --- Branch Reverb        
        BranchReverb->processAudio(&pitchOut, bran_rev_out);

        // --- Master Reverb        
        // Mix branch reverb output with dry input
        mast_rev_in[0] = shim_shimmer * bran_rev_out[0] + (1 - shim_shimmer) * inL[i];
        mast_rev_in[1] = shim_shimmer * bran_rev_out[1] + (1 - shim_shimmer) * inR[i];

        // Process master reverb        
        MasterReverb->processAudio(mast_rev_in, mast_rev_out);        

        // Stereo spread processing + output allocation
        outL[i] = _wet * mast_rev_out[0] + _dry * inL[i];
        outR[i] = _wet * mast_rev_out[1] + _dry * inR[i];
    }

   // Write samples to file
   //string post = "test_output.txt";
   //WriteBufferToFile(outputs, sampleFrames, post);

}
/*--------------------------------------------------------------------*/



/* ------------------------------------------------------------------------------------------------------------
  ------------------------------------------  PARAMETERS  ------------------------------------------------------
  ------------------------------------------------------------------------------------------------------------ */
  // set reverb parameters values
void Shimmer::setParameter(VstInt32 index, float value)
{
    switch (index) {
    case Param_mix: {
        shim_mix = value;        
        updateMix();
        break;
    }
    case Param_roomSize: {
        shim_roomSize = value;
        BranchReverb->setRoomSize(shim_roomSize);
        MasterReverb->setRoomSize(shim_roomSize);
        break;
    }
    case Param_shimmer: {
        shim_shimmer = value;
        break;
    }
    case Param_decay: {
        shim_decay = value;
        BranchReverb->setDecayInSeconds(0.25 * shim_decay * MAX_REVERB_DECAY_IN_SECONDS);
        MasterReverb->setDecayInSeconds(shim_decay * MAX_REVERB_DECAY_IN_SECONDS);
        break;
    }    
    case Param_damping: {
        shim_damping = value; 
        float freq = exp(mapValueIntoRange(1.0 - shim_damping, MIN_DAMPING_FREQUENCY_LOG, MAX_DAMPING_FREQUENCY_LOG));
        BranchReverb->setDampingFrequency(freq);
        MasterReverb->setDampingFrequency(freq);
        break;
    }        
    case Param_space: {
        shim_space = value;
        MasterReverb->setStereoSpread(shim_space);
        break;
    }
    case Param_intervals: {        
        if (value == 1)
            value = 0.99; // if value = 1, then pitIdx = NUM_OF_PITCH_INTRVL_ALLOWED + 1 -> outside of array boundaries
        shim_intervals = value;
        int pitIdx = shim_intervals / DELTA_PARAMETER_BETWEEN_INTERVALS;
        PitchShift_1oct->setPitchShift(INTERVALS_IN_SEMITONES_PITCH1[pitIdx]);
        PitchShift_2oct->setPitchShift(INTERVALS_IN_SEMITONES_PITCH2[pitIdx]);
        updateMixPitchShifters(INTERVALS_IN_SEMITONES_PITCH2[pitIdx]);
        break;
    }
    case Param_modDepth: {
        shim_modDepth = value;
        MasterReverb->setModDepth(shim_modDepth);
        break;
    }
    case Param_modRate: {
        shim_modRate = value;
        MasterReverb->setModRate(shim_modRate * MAX_MOD_RATE);
        break;
    }
    case Param_lpf: {
        shim_lpf = value;
        float lpf = exp(mapValueIntoRange(value, LPF_FILTER_MIN_FREQ_LOG, LPF_FILTER_MAX_FREQ_LOG));
        MasterReverb->setLowPassFrequency(lpf);
        break;
    }
    case Param_hpf: {        
        shim_hpf = value;
        float hpf = mapValueIntoRange(value, HPF_FILTER_MIN_FREQ, HPF_FILTER_MAX_FREQ);
        MasterReverb->setHighPassFrequency(hpf);
        break;
    }    
    default:
        break;
    }
}
/*--------------------------------------------------------------------*/

/*--------------------------------------------------------------------*/
// return reverb parameters values
float Shimmer::getParameter(VstInt32 index)
{
    float param = 0;
    switch (index) {
    case Param_mix:
    {
        param = shim_mix;
        break;
    }
    case Param_roomSize: 
    {
        param = shim_roomSize;
        break;
    }
    case Param_shimmer:
    {
        param = shim_shimmer;
        break;
    }
    case Param_decay:
    {
        param = shim_decay;
        break;
    }    
    case Param_damping:
    {
        param = shim_damping;
        break;
    }    
    case Param_space:
    {
        param = shim_space;
        break;
    }
    case Param_intervals:
    {
        param = shim_intervals;
        break;
    }
    case Param_modDepth: 
    {
        param = shim_modDepth;
        break;
    }
    case Param_modRate: 
    {
        param = shim_modRate;
        break;
    }
    case Param_lpf: 
    {
        //param = mapValueOutsideRange(log(shim_lpf), MIN_LPF_FREQUENCY_LOG, MAX_LPF_FREQUENCY_LOG);
        param = shim_lpf;
        break;
    }
    case Param_hpf: 
    {        
        //param = mapValueOutsideRange(shim_hpf, HPF_FILTER_MIN_FREQ, HPF_FILTER_MAX_FREQ);
        param = shim_hpf;
        break;
    }    
    default:
        break;
    }
    return param;
}
/*--------------------------------------------------------------------*/

/*--------------------------------------------------------------------*/
// return reverb parameters labels
void Shimmer::getParameterLabel(VstInt32 index, char* label)
{
    switch (index) {
    case Param_mix: {
        vst_strncpy(label, "", kVstMaxParamStrLen);
        break;
    }
    case Param_roomSize: {
        vst_strncpy(label, "", kVstMaxParamStrLen);
        break;
    }
    case Param_shimmer: {
        vst_strncpy(label, "", kVstMaxParamStrLen);
        break;
    }
    case Param_decay: {
        vst_strncpy(label, "s", kVstMaxParamStrLen);
        break;
    }
    case Param_damping: {
        vst_strncpy(label, "", kVstMaxParamStrLen);
        break;
    }
    case Param_space: {
        vst_strncpy(label, "", kVstMaxParamStrLen);
        break;
    }
    case Param_intervals: {
        vst_strncpy(label, "", kVstMaxParamStrLen);
        break;
    }
    case Param_modRate: {
        vst_strncpy(label, "Hz", kVstMaxParamStrLen);
        break;
    }
    case Param_modDepth: {
        vst_strncpy(label, "", kVstMaxParamStrLen);
        break;
    }
    case Param_lpf: {
        vst_strncpy(label, "Hz", kVstMaxParamStrLen);
        break;
    }
    case Param_hpf: {
        vst_strncpy(label, "Hz", kVstMaxParamStrLen);
        break;
    }      
    default: {
        break;
    }
    }
}
/*--------------------------------------------------------------------*/

/*--------------------------------------------------------------------*/
// return reverb parameters values for displaying purposes
void Shimmer::getParameterDisplay(VstInt32 index, char* text)
{    
    switch (index) {
    case Param_mix: {
        float2string(shim_mix * 10, text, kVstMaxParamStrLen);
        break;
    }
    case Param_roomSize: {
        float2string(shim_roomSize * 10, text, kVstMaxParamStrLen);
        break;
    }
    case Param_shimmer: {
        float2string(shim_shimmer * 10, text, kVstMaxParamStrLen);
        break;
    }
    case Param_decay: {
        float2string(shim_decay * MAX_REVERB_DECAY_IN_SECONDS, text, kVstMaxParamStrLen);
        break;
    }
    case Param_damping: {
        float2string(shim_damping * 10, text, kVstMaxParamStrLen);
        break;
    }
    case Param_space: {
        float2string(shim_space * 10, text, kVstMaxParamStrLen);
        break;
    }
    case Param_intervals: {
        int pitIdx = shim_intervals / DELTA_PARAMETER_BETWEEN_INTERVALS;
        vst_strncpy(text, INTERVALS_NAMES_STRING[pitIdx], kVstMaxParamStrLen);
        break; 
    }
    case Param_modRate: {
        float2string(shim_modRate * MAX_MOD_RATE, text, kVstMaxParamStrLen);
        break;
    }
    case Param_modDepth: {
        float2string(shim_modDepth * 10, text, kVstMaxParamStrLen);
        break;
    }
    case Param_lpf: {
        float lpf = exp(mapValueIntoRange(shim_lpf, LPF_FILTER_MIN_FREQ_LOG, LPF_FILTER_MAX_FREQ_LOG));
        float2string(lpf, text, kVstMaxParamStrLen);
        break;
    }
    case Param_hpf: {
        float hpf = mapValueIntoRange(shim_hpf, HPF_FILTER_MIN_FREQ, HPF_FILTER_MAX_FREQ);
        float2string(hpf, text, kVstMaxParamStrLen);
        break;
    }      
    default: {
        break;
    }
    }
}
/*--------------------------------------------------------------------*/

/*--------------------------------------------------------------------*/
// return the parameter's name for displaying purpose
void Shimmer::getParameterName(VstInt32 index, char* text)
{
    switch (index) {
    case Param_mix: {
        vst_strncpy(text, "Mix", kVstMaxParamStrLen);
        break;
    }
    case Param_roomSize: {
        vst_strncpy(text, "Size", kVstMaxParamStrLen);
        break;
    }
    case Param_shimmer: {
        vst_strncpy(text, "Shimmer", kVstMaxParamStrLen);
        break;
    }
    case Param_decay: {
        vst_strncpy(text, "Decay", kVstMaxParamStrLen);
        break;
    }
    case Param_damping: {
        vst_strncpy(text, "Damping", kVstMaxParamStrLen);
        break;
    }
    case Param_space: {
        vst_strncpy(text, "Space", kVstMaxParamStrLen);
        break;
    }
    case Param_intervals: {
        vst_strncpy(text, "Intervals", kVstMaxParamStrLen);
        break;
    }
    case Param_modRate: {
        vst_strncpy(text, "Mod Rate", kVstMaxParamStrLen);
        break;
    }
    case Param_modDepth: {
        vst_strncpy(text, "Mod Depth", kVstMaxParamStrLen);
        break;
    }
    case Param_lpf: {
        vst_strncpy(text, "LPF", kVstMaxParamStrLen);
        break;
    }
    case Param_hpf: {
        vst_strncpy(text, "HPF", kVstMaxParamStrLen);
        break;
    }     
    default: {
        break;
    }
    }
}
/*--------------------------------------------------------------------*/


/* ------------------------------------------------------------------------------------------------------------
 -------------------------------------------  PROGRAM  --------------------------------------------------------
 ------------------------------------------------------------------------------------------------------------ */

 /*--------------------------------------------------------------------*/
 // Define presets parameters values
void Shimmer::InitPresets()
{
    // Instruction to allow the DAW to save parameters as currently shown on the UI
    programsAreChunks(true);

    // Create object for presets
    shim_presets = new ShimmerPresets[NUM_PRESETS];

    /*----------------------------------------------------*/
    // "Default" preset
    strcpy(shim_presets[0].name, "Default");    
    shim_presets[0].shim_mix = 0.5;
    shim_presets[0].shim_roomSize = 0.4;
    shim_presets[0].shim_decay = 0.2;
    shim_presets[0].shim_damping = 0.5;
    shim_presets[0].shim_space = 0.5;
    shim_presets[0].shim_shimmer = 0.5;
    shim_presets[0].shim_intervals = 0.71;
    shim_presets[0].shim_modRate = 0.0;
    shim_presets[0].shim_modDepth = 0.0;
    shim_presets[0].shim_lpf = 1.0;
    shim_presets[0].shim_hpf = 0.0;

    /*----------------------------------------------------*/
    // "Small Room" 
    strcpy(shim_presets[1].name, "Small Room");
    shim_presets[1].shim_mix = 0.5;
    shim_presets[1].shim_roomSize = 0.1;
    shim_presets[1].shim_decay = 0.1;
    shim_presets[1].shim_damping = 0.15;
    shim_presets[1].shim_space = 0.25;
    shim_presets[1].shim_shimmer = 0.0;
    shim_presets[1].shim_intervals = 0.1;
    shim_presets[1].shim_modRate = 0.0;
    shim_presets[1].shim_modDepth = 0.0;
    shim_presets[1].shim_lpf = 1.0;
    shim_presets[1].shim_hpf = 0.0;
     
    // "Hall" 
    strcpy(shim_presets[2].name, "Hall");
    shim_presets[2].shim_mix = 0.7;
    shim_presets[2].shim_roomSize = 0.6;
    shim_presets[2].shim_decay = 0.4;
    shim_presets[2].shim_damping = 0.3;
    shim_presets[2].shim_space = 0.6;
    shim_presets[2].shim_shimmer = 0.0;
    shim_presets[2].shim_intervals = 0.1;
    shim_presets[2].shim_modRate = 0.0;
    shim_presets[2].shim_modDepth = 0.0;
    shim_presets[2].shim_lpf = mapValueOutsideRange(log(16000.0), MIN_LPF_FREQUENCY_LOG, MAX_LPF_FREQUENCY_LOG);
    shim_presets[2].shim_hpf = 0.0;

    // "Ambience Damped" 
    strcpy(shim_presets[3].name, "Ambience Damped");
    shim_presets[3].shim_mix = 0.8;
    shim_presets[3].shim_roomSize = 0.6;
    shim_presets[3].shim_decay = 0.6;
    shim_presets[3].shim_damping = 0.8;
    shim_presets[3].shim_space = 0.7;
    shim_presets[3].shim_shimmer = 1.0;
    shim_presets[3].shim_intervals = 0.71;
    shim_presets[3].shim_modRate = 0.0;
    shim_presets[3].shim_modDepth = 0.0;
    shim_presets[3].shim_lpf = 1.0;
    shim_presets[3].shim_hpf = mapValueOutsideRange(150.0, MIN_HPF_FREQUENCY, MAX_HPF_FREQUENCY);

    // "Ambience Modulated" 
    strcpy(shim_presets[4].name, "Ambience Modulated");
    shim_presets[4].shim_mix = 0.8;
    shim_presets[4].shim_roomSize = 0.6;
    shim_presets[4].shim_decay = 0.7;
    shim_presets[4].shim_damping = 0.2;
    shim_presets[4].shim_space = 0.7;
    shim_presets[4].shim_shimmer = 1.0;
    shim_presets[4].shim_intervals = 0.81;
    shim_presets[4].shim_modRate = 0.4;
    shim_presets[4].shim_modDepth = 0.5;
    shim_presets[4].shim_lpf = mapValueOutsideRange(log(17000.0), MIN_LPF_FREQUENCY_LOG, MAX_LPF_FREQUENCY_LOG);;
    shim_presets[4].shim_hpf = 0.0;    

    // Set the program when creating a new plugin instance
    int initIdx = 0;
    Shimmer::setProgram(initIdx);
    shim_mix = shim_presets[initIdx].shim_mix;
    shim_roomSize = shim_presets[initIdx].shim_roomSize;
    shim_decay = shim_presets[initIdx].shim_decay;
    shim_damping = shim_presets[initIdx].shim_damping;
    shim_space = shim_presets[initIdx].shim_space;
    shim_shimmer = shim_presets[initIdx].shim_shimmer;
    shim_intervals = shim_presets[initIdx].shim_intervals;
    shim_modRate = shim_presets[initIdx].shim_modRate;
    shim_modDepth = shim_presets[initIdx].shim_modDepth;
    shim_lpf = shim_presets[initIdx].shim_lpf;
    shim_hpf = shim_presets[initIdx].shim_hpf;
}
/*--------------------------------------------------------------------*/

/*--------------------------------------------------------------------*/
void Shimmer::setProgram(VstInt32 program)
{
    // Call the current implementation of "setProgram"
    AudioEffect::setProgram(program);

    // Create an instante of ShimmerPresets with current preset
    ShimmerPresets* cp = &shim_presets[program];

    // Set each parameter
    setParameter(Param_mix, cp->shim_mix);
    setParameter(Param_roomSize, cp->shim_roomSize);
    setParameter(Param_decay, cp->shim_decay);
    setParameter(Param_damping, cp->shim_damping);
    setParameter(Param_space, cp->shim_space);
    setParameter(Param_shimmer, cp->shim_shimmer);
    setParameter(Param_intervals, cp->shim_intervals);
    setParameter(Param_modRate, cp->shim_modRate);
    setParameter(Param_modDepth, cp->shim_modDepth);
    setParameter(Param_lpf, cp->shim_lpf);
    setParameter(Param_hpf, cp->shim_hpf);
}
/*--------------------------------------------------------------------*/

/*--------------------------------------------------------------------*/
// return program's current name
void Shimmer::getProgramName(char* name)
{
    strcpy(name, shim_presets[curProgram].name);
}
/*--------------------------------------------------------------------*/

/*--------------------------------------------------------------------*/
// get the program's list
bool Shimmer::getProgramNameIndexed(VstInt32 category, VstInt32 index, char* text)
{
    if (index < NUM_PRESETS)
    {
        strcpy(text, shim_presets[index].name);
        return true;
    }
    return false;
}
/*--------------------------------------------------------------------*/

/* ------------------------------------------------------------------------------------------------------------
 ------------------------------------------  SAVE&LOAD  -------------------------------------------------------
 ------------------------------------------------------------------------------------------------------------ */

// Save current parameters
VstInt32 Shimmer::getChunk(void** data, bool isPreset) 
{    
    ShimmerParameters params;
    params.mix = this->shim_mix;
    params.size = this->shim_roomSize;
    params.decay = shim_decay;
    params.shimmer = shim_shimmer;
    params.interval = shim_intervals;
    params.damping = shim_damping;
    params.rate = shim_modRate;
    params.depth = shim_modDepth;
    params.space = shim_space;
    params.lpf = shim_lpf;
    params.hpf = shim_hpf;   

    *data = &params;

    return sizeof(params);    
}

// Load parameters saved
VstInt32 Shimmer::setChunk(void* data, VstInt32 byteSize, bool isPreset)
{
    ShimmerParameters* cp = (ShimmerParameters*)data;
    setParameter(Param_mix, cp->mix);
    setParameter(Param_roomSize, cp->size);
    setParameter(Param_decay, cp->decay);
    setParameter(Param_damping, cp->damping);
    setParameter(Param_space, cp->space);
    setParameter(Param_shimmer, cp->shimmer);
    setParameter(Param_intervals, cp->interval);
    setParameter(Param_modRate, cp->rate);
    setParameter(Param_modDepth, cp->depth);
    setParameter(Param_lpf, cp->lpf);
    setParameter(Param_hpf, cp->hpf);

    return 0;
}

/* ------------------------------------------------------------------------------------------------------------
 ---------------------------------------------  NAME  ---------------------------------------------------------
 ------------------------------------------------------------------------------------------------------------ */

 /*--------------------------------------------------------------------*/
bool Shimmer::getEffectName(char* name)
{
    vst_strncpy(name, "Shimmer", kVstMaxEffectNameLen);
    return true;
}
/*--------------------------------------------------------------------*/

/*--------------------------------------------------------------------*/
bool Shimmer::getVendorString(char* name)
{
    vst_strncpy(name, "Fox Suite", kVstMaxVendorStrLen);
    return true;
}
/*--------------------------------------------------------------------*/

/* ------------------------------------------------------------------------------------------------------------
 ------------------------------------------  DESTRUCTOR  ------------------------------------------------------
 ------------------------------------------------------------------------------------------------------------ */
Shimmer::~Shimmer()
{
    //Free BranchReverb, MasterReverb and pitch shifters
    delete MasterReverb;
    delete BranchReverb;
    delete PitchShift_1oct, PitchShift_2oct;
}


