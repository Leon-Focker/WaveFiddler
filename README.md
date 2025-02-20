# WaveFiddler
WaveFiddler is a command line application for creative conversions between audio and image files.

Precompiled binaries can be found in the [Releases tab](https://github.com/Leon-Focker/WaveFiddler/releases/)

## Currently supported conversions:
- From .wav to one image (.jpg or .png) with variable FFT size.
- From .wav to an image sequence, which will be in sync with the audio file when played at the given frame rate (can be set).
- From an image (.jpg or .png) to a wavetable (.wav) with one frame.
- From an image (.jpg or .png) to a wavetable (.wav) with multiple frames.

## Examples
To run these examples, make sure to adjust the paths to the WaveFiddler binary and this examples folder based on their locations on your system. Run ``./wavefiddler --help`` for detailed documentation on usage and available options. 

**Sound to Image**

Generating a black and white image sequence from a mono file.

``.\wavefiddler.exe .\examples\kick.wav -a 1``

Changing the color map to get a different hue, only generate the first frame

``.\wavefiddler.exe .\examples\kick.wav -a 1 -s 0 -c 0.2,0,1,0,0.5,0``

Use a different method to generate a single image

``.\wavefiddler.exe .\examples\bells.wav -s 0``

Use yet another method, start at the 24000s sample of the soundfile, use 40000 samples and make it vibrant

``.\wavefiddler.exe .\examples\bells.wav -s 24000 -a 4 -f 40000 -v 5 -c 0.5,0,0,1,0,0.2``

**Image to sound**

Get a wavetable of length 512 from an image, plot the generated waveform and spectrum

``.\wavefiddler.exe .\examples\cat.jpg --plot-waveforms --plot-spectra``

Get a wavetable of length 512 from an image, plot the generated waveform and spectrum but fit the original spectrum into the length of the wavetable

``.\wavefiddler.exe .\examples\cat.jpg -S --plot-waveforms --plot-spectra``

Get a wavetable with 100 frames and automatic fft-size, name the output file

``.\wavefiddler.exe .\examples\cat.jpg -f 0 -i 100 -n "cat.wav"``

``.\wavefiddler.exe .\examples\phillippa.png -f 0 -i 100 -n "phillippa"``