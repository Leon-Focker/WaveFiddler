# WaveFiddler
WaveFiddler is a command line application for creative conversions between audio and image files.

Precompiled binaries can be found in the [Releases tab](https://github.com/Leon-Focker/WaveFiddler/releases/)

Since this software is not signed for MacOS, you might have to do something like this:<br>
``sudo chmod a+x ./wavefiddler``<br>
MacOS < 15.0: ``--add ./wavefiddler`` <br>
MacOS >= 15.0 ``-d com.apple.quarantine ./wavefiddler``

## Currently supported conversions:
- From .wav to one image (.jpg or .png) with variable FFT size.
- From .wav to an image sequence, which will be in sync with the audio file when played at the given frame rate (can be set).
- From an image (.jpg or .png) to a wavetable (.wav) with one frame.
- From an image (.jpg or .png) to a wavetable (.wav) with multiple frames.

## Examples
To run these examples, make sure to adjust the paths to the WaveFiddler binary and this examples folder based on their locations on your system. When you don't specify an output directory, the results are put into ./wavefiddler_outputs/ (which will be created if it does not exist). Run ``./wavefiddler --help`` for detailed documentation on usage and available options. 

**Sound to Image**

Generating a black and white image sequence from a mono file.<br>
``.\wavefiddler.exe .\examples\kick.wav -a 1``

https://github.com/user-attachments/assets/7773feac-a899-4f7e-a0a1-66734812a3fa

Changing the color map to get a different hue, only generate the first frame<br>
``.\wavefiddler.exe .\examples\kick.wav -a 1 -s 0 -c 0.2,0,1,0,0.5,0``

![kick](https://github.com/user-attachments/assets/3765e391-443a-43d4-934c-928efd2dc602)

Specifying image dimensions and using a different method<br>
``.\wavefiddler.exe .\examples\kick.wav -a 3 -s 0 -c 0.2,0,1,0,0.5,0 -d 400,200``

![kicktt](https://github.com/user-attachments/assets/9096561e-52c1-464c-b6b6-f87eb9c9e15a)

Use a different method to generate a single image. A stereo soundfile produces more interesting colours<br>
``.\wavefiddler.exe .\examples\bells.wav -s 0``

![bells1](https://github.com/user-attachments/assets/ffe9e2bf-3c3d-4fe4-a774-cacfa56aa64c)

Use yet another method, start at the 24000s sample of the soundfile, use 40000 samples and make it vibrant<br>
``.\wavefiddler.exe .\examples\bells.wav -s 24000 -a 4 -f 40000 -v 5 -c 0.5,0,0,1,0,0.2``

![bells](https://github.com/user-attachments/assets/bd739e06-2808-4dac-bec2-a0b926e4e47c)

**Image to sound**

Get a wavetable of length 512 from an image, plot the generated waveform and spectrum<br>
``.\wavefiddler.exe .\examples\cat.jpg --plot-waveforms --plot-spectra``

![example_12](https://github.com/user-attachments/assets/8b380739-405e-4079-a3d8-8f730c753955)
![example_11](https://github.com/user-attachments/assets/ef0a774f-8c28-4ae2-82bd-47261d160847)

Get a wavetable of length 512 from an image, plot the generated waveform and spectrum but fit the original spectrum into the length of the wavetable<br>
``.\wavefiddler.exe .\examples\cat.jpg -S --plot-waveforms --plot-spectra``

![waveform](https://github.com/user-attachments/assets/013add22-ce85-48f3-856a-83da052e2299)
![spectrum_grey](https://github.com/user-attachments/assets/2748740c-0a77-4766-b7a2-f8b032daff63)

Get a wavetable with 100 frames and automatic fft-size, name the output file<br>
``.\wavefiddler.exe .\examples\cat.jpg -f 0 -i 100 -n "cat.wav"``

https://github.com/user-attachments/assets/7298ed7b-38be-4c8e-8a08-6db47ced0aca

``.\wavefiddler.exe .\examples\phillippa.png -f 0 -i 100 -n "phillippa"``

https://github.com/user-attachments/assets/b38eb19e-9ad5-4f1d-9d55-a7cf444717d6
