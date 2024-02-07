import("stdfaust.lib");

waveformDisplay(index, oscillator, freq) = hgroup("[%index]Wave %index", oscillator(freq*multiplier)*gain_wf)
with{
    gain_wf = hslider("[0]gain %index[style:knob]",1, 0, 1, 0.01);
    multiplier = hslider("[1] frequency multiplier %index", 1 , 0, 4, 0.01);
};
//ovo je u sustini za gui, za svaku od komponenata signala, ovde konkretno je dostupan gain i multiplikator frekvencije


subtractive(index, oscillator, freq) = vgroup("[%index]Waveform %index", waveformDisplay(index, oscillator, freq) : hgroup("[1]Filter", fi.resonlp(resFreq, q ,1))
with{
    ctFreq = hslider("[0]Cutoff Frequency[style:knob]", 2000,50,10000,0.1);
    q = hslider("[1]Q[style:knob]",5,1,30,0.1);
    lfoFreq = hslider("[2]LFO Frequency[style:knob]",10,0.1,20,0.01);
    lfoDepth = hslider("[3]LFO Depth[style:knob]",500,1,10000,1);
    resFreq = os.osc(lfoFreq)*lfoDepth + ctFreq: max(30);
    
});
//svaka od komponenata ima primenjen resonant low pass filter


envelope = hgroup("[5]Envelope", en.adsr(attack,decay,sustain,release,gate)*gain*0.1)
with{
    attack = hslider("[0]Attack[style:knob]",50,1,1000,1)*0.001;
    decay = hslider("[1]Decay[style:knob]",50,1,1000,1)*0.001;
    sustain = hslider("[2]Sustain[style:knob]",0.8,0.01,1,0.001)*0.001;
    release = hslider("[3]Release[style:knob]",50,1,1000,1)*0.001;
    gain = hslider("[4]Gain[style:knob]",1,0,1,0.01);
    gate = button("[5]gate");
};
//ADSR (Attack, Decay, Sustain, Release) envelope generator.


noise(freq) = no.noise;
signal(freq) = (subtractive(0,noise,freq) + subtractive(1,os.triangle,freq) + subtractive(2, os.square, freq) + subtractive(3, os.sawtooth, freq)) / 4 ;
//krajnji signal se sastoji od suma i 3 oscilatora


reverb = on,_,_,hgroup("Reverb", re.stereo_freeverb(fb1 ,fb2, damp, spread)) :> _,_,_,_,_
with{
    on = nentry("Add Reverb", 1, 0, 1, 1);
    fb1 = hslider ("[1]fb1[style:knob]", 0.25, 0 , 1, 0.01);
    fb2 = hslider("[2]fb2[style:knob]",0.25, 0 , 1, 0.01);
    damp = hslider("[3]damp[style:knob]",0.25, 0 , 1, 0.01);
    spread = hslider("[4]spread", 4, 1, 16, 1);
};
//dodat audio efekat reverb
//A simple Schroeder reverberator primarily developed by "Jezar at Dreampoint" that is extensively used in the free-software world. It uses four Schroeder allpasses in series and eight parallel Schroeder-Moorer filtered-feedback comb-filters for each audio channel, and is said to be especially well tuned.

add_reverb(rev, ch1, ch2, ch3, ch4) = select2(rev, ch1, ch3),select2(rev,ch2,ch4);
//2-4 mapiranje kod reverba, valjda

process = vgroup("Subtractive Synthesizer", signal(freq)*envelope) <: (reverb: add_reverb) 
with {
    freq= hslider("[7]freq[style:knob]",440, 50, 2000, 0.01);
};