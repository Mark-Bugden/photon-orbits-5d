// The trail and the photon. The trail is a thick, unlit line, so its colour carries only ψ, with a thin dark edge
// that keeps it legible where it crosses the spheroid's outline. Its buffers are allocated once; each frame,
// trail.js's fill() rewrites the part the trail uses, oldest step first.
//
// Every strand lies on the same surface, so where two cross they are at the same depth, and a depth test picks a
// winner pixel by pixel: a jagged join, with no edge between them. So the trail doesn't test depth against itself.
// It is drawn twice, once for each face of the spheroid: the steps on the far face, then those on the near face,
// each pass in time order. Each step draws its edge and its colour together, so the newest strand passes cleanly
// over older ones, edge and all, and the near face always covers the far one. The far pass still tests depth, so
// the horizon hides the strands behind it, and the near pass writes depth, so the translucent shell, drawn last,
// dims only the far strands.
//
// The lines are smoothed by multisampling alone, without alphaToCoverage: three's canvas always keeps an alpha
// channel, and the soft rims that alphaToCoverage gives the round ends are written as partial transparency, which
// let a light page show through as faint rings at every joint.

import {
  AlwaysDepth, BackSide, Color, DynamicDrawUsage, Group, InstancedInterleavedBuffer, InterleavedBufferAttribute, Mesh,
  MeshBasicMaterial, SphereGeometry, Vector3,
} from 'three';
import { LineMaterial } from 'three/addons/lines/LineMaterial.js';
import { LineSegments2 } from 'three/addons/lines/LineSegments2.js';
import { LineSegmentsGeometry } from 'three/addons/lines/LineSegmentsGeometry.js';

import { WHEEL, wheelIndex } from '../colour.js';
import { fill } from '../trail.js';
import { BEHIND } from './surfaces.js';

export const CAPACITY = 65536; // segments: a trail of λ = 300 needs about 20,000, or more near the extremal edge
const WIDTH = 3; // CSS pixels of colour
const EDGE = 1; // CSS pixels of dark edge each side
const FADED = 0.15; // the colour left at the tail, as a fraction of the way from BEHIND to the full colour
const FULL = 0.5; // the newest half of the trail is at full colour

// An interleaved, instanced buffer of 6 floats a segment (start, then end), and its two attributes.
function segments(geometry, array, start, end) {
  const buffer = new InstancedInterleavedBuffer(array, 6, 1).setUsage(DynamicDrawUsage);
  geometry.setAttribute(start, new InterleavedBufferAttribute(buffer, 3, 0));
  geometry.setAttribute(end, new InterleavedBufferAttribute(buffer, 3, 3));
  return buffer;
}

// Writes the colour of the angle ψ + shift, faded towards BEHIND by its place along the trail, in linear RGB.
const behind = new Color(BEHIND);
let shift = 0;
function paint(psi, fraction, out, offset) {
  const k = 3 * wheelIndex(psi + shift);
  const weight = FADED + (1 - FADED) * Math.min(1, fraction / FULL);
  out[offset] = behind.r + weight * (WHEEL[k] - behind.r);
  out[offset + 1] = behind.g + weight * (WHEEL[k + 1] - behind.g);
  out[offset + 2] = behind.b + weight * (WHEEL[k + 2] - behind.b);
}

// Changes to three's line shader (r187): draw only the steps on one face of the spheroid, and paint the outer part
// of the line's width as the edge. A step is on the near face when the spheroid's outward normal at its middle,
// (X/A², Y/A², Z/B²), points towards the camera. Distance across the line is |vUv.x|, which in the round end caps
// is the distance across the step rather than from its end, so a step's cap matches the edge of the step it
// overlaps. Each change throws if three's shader no longer has the line it hooks onto.
const PATCHES = {
  vertexShader: [
    ['attribute vec3 instanceColorEnd;', 'attribute vec3 instanceColorEnd;\nuniform vec3 axes;\nuniform float side;'],
    ['gl_Position = clip;', `gl_Position = clip;
      vec3 middle = 0.5 * ( instanceStart + instanceEnd );
      if ( side * dot( middle / axes, cameraPosition - middle ) < 0.0 ) gl_Position = vec4( 0.0, 0.0, 2.0, 1.0 ); // off this face: outside the view, so not drawn`],
  ],
  fragmentShader: [
    ['uniform float linewidth;', 'uniform float linewidth;\nuniform vec3 edge;\nuniform float inner;'],
    ['#include <color_fragment>', `#include <color_fragment>
      float across = abs( vUv.x );
      float blur = 0.5 * fwidth( across ); // a one-pixel blend between the colour and the edge
      diffuseColor.rgb = mix( diffuseColor.rgb, edge, smoothstep( inner - blur, inner + blur, across ) );`],
  ],
};

function patch(shader) {
  for (const [stage, changes] of Object.entries(PATCHES)) {
    for (const [anchor, replacement] of changes) {
      if (!shader[stage].includes(anchor)) throw new Error(`three's line shader has changed: no "${anchor}"`);
      shader[stage] = shader[stage].replace(anchor, replacement);
    }
  }
}

/** { group, setSpheroid(axes), show(path, head, length, psi0), hide(), photon }: the trail and the photon. */
export function createPath() {
  const positions = new Float32Array(6 * CAPACITY);
  const colours = new Float32Array(6 * CAPACITY);
  const geometry = new LineSegmentsGeometry();
  const positionBuffer = segments(geometry, positions, 'instanceStart', 'instanceEnd');
  const colourBuffer = segments(geometry, colours, 'instanceColorStart', 'instanceColorEnd');
  geometry.instanceCount = 0;

  const axes = new Vector3(1, 1, 1); // the spheroid's squared semi-axes, shared by both passes
  const face = (side, depth) => {
    const material = new LineMaterial({ vertexColors: true, linewidth: WIDTH + 2 * EDGE, fog: true, ...depth });
    Object.assign(material.uniforms, {
      axes: { value: axes },
      side: { value: side },
      edge: { value: new Color(0x07090c) },
      inner: { value: WIDTH / (WIDTH + 2 * EDGE) },
    });
    material.onBeforeCompile = patch;
    const mesh = new LineSegments2(geometry, material);
    mesh.frustumCulled = false; // the buffers change every frame
    return mesh;
  };
  const far = face(-1, { depthWrite: false });
  const near = face(1, { depthFunc: AlwaysDepth });
  far.renderOrder = 2;
  near.renderOrder = 3;

  // The photon: a white dot, ringed in black by the inside of a slightly larger sphere, drawn after the trail.
  const ball = new SphereGeometry(1, 24, 16);
  const ringed = new Mesh(ball, new MeshBasicMaterial({ color: 0x000000, side: BackSide }));
  ringed.scale.setScalar(1.5);
  const photon = new Group().add(new Mesh(ball, new MeshBasicMaterial({ color: 0xffffff })), ringed);
  photon.renderOrder = 4;

  const group = new Group().add(far, near, photon);
  const head = [0, 0, 0, 0];

  return {
    group,
    photon,
    /** The orbit's spheroid, as its semi-axes [across, up], which tells the two passes which face a step is on. */
    setSpheroid([across, up]) {
      axes.set(across * across, across * across, up * up);
    },
    /**
     * Draws the trail from λ = head − length to the photon at λ = head, coloured by ψ + psi0, where psi0 is ψ at
     * the start of the path. Returns the photon's ψ + psi0.
     */
    show(path, lam, length, psi0 = 0) {
      shift = psi0;
      const count = fill(path, lam, length, positions, colours, paint, head);
      geometry.instanceCount = count;
      for (const buffer of [positionBuffer, colourBuffer]) {
        buffer.clearUpdateRanges();
        buffer.addUpdateRange(0, 6 * count);
        buffer.needsUpdate = true;
      }
      photon.position.set(head[0], head[1], head[2]);
      group.visible = true;
      return head[3] + psi0;
    },
    hide() {
      group.visible = false;
    },
  };
}
