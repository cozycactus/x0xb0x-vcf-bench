# x0xb0x VCF Bench

Play it: **https://cozycactus.github.io/x0xb0x-vcf-bench/**

A real-time model of the x0xb0x filter (the TB-303-style four-pole diode ladder), with its envelope, Env Mod, Accent and Decay circuits, the VCA, the saw and square oscillator, and a 16-step sequencer with accent and slide. It runs in the browser as an AudioWorklet and works on a phone.

The model copies an ngspice simulation of one x0xb0x build's KiCad schematic, which in turn copies the circuit. It is a model of a model: it was fitted to that simulation, not to a real TB-303 or x0xb0x, and no real unit has been measured. Against the simulation, it lands within 0.027 octave on the filter's peak frequency, within 0.13 dB RMS (median over 25 knob settings) on the response, and within 0.036 octave on the envelope and accent peaks.

## Playing

Press Play, then drag the knobs up and down (arrow keys work too; double-click resets a knob). Tap a step and then a key to set its note; Down, Up, Rest, Accent and Slide change the selected step. If a slow phone crackles, switch "4x" to "2x".

The Space panel is a reverb added on the page (Space, Size, Tone). The x0xb0x has no reverb; the model of the circuit is `dsp.js`, and the reverb sits outside it.

## Running it locally

The page loads `dsp.js` as an AudioWorklet module, so it needs to be served rather than opened as a file:

```bash
python3 -m http.server 8000
```

then open http://localhost:8000.

## Files

- `index.html`: the page.
- `dsp.js`: the model and the AudioWorklet processor.

The ngspice study, the fitting and the verification scripts live in a separate, private repository. `dsp.js` is CC0.
