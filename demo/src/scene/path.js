// The trail and the photon. The trail is a thick, unlit line, so its colour carries only ψ, in a thin dark casing
// that keeps it legible where it crosses the spheroid's outline. Its buffers are allocated once; each frame,
// trail.js's fill() rewrites the part the trail uses.
//
// The line is drawn first and marks its pixels in the stencil buffer, and the casing is drawn only where it
// isn't. A depth offset can't do this: each step's round end caps are flat in depth, so where a step heads away
// from the camera, the next step's casing cap sits in front of the line it overlaps, and the line looks dashed.
//
// The lines are smoothed by multisampling alone, without alphaToCoverage: three's canvas always keeps an alpha
// channel, and the soft rims that alphaToCoverage gives the round ends are written as partial transparency, which
// let a light page show through as faint rings at every joint.

import {
  AlwaysStencilFunc, BackSide, Color, DynamicDrawUsage, Group, InstancedInterleavedBuffer, InterleavedBufferAttribute,
  Mesh, MeshBasicMaterial, NotEqualStencilFunc, ReplaceStencilOp, SphereGeometry,
} from 'three';
import { LineMaterial } from 'three/addons/lines/LineMaterial.js';
import { LineSegments2 } from 'three/addons/lines/LineSegments2.js';
import { LineSegmentsGeometry } from 'three/addons/lines/LineSegmentsGeometry.js';

import { WHEEL, wheelIndex } from '../colour.js';
import { fill } from '../trail.js';
import { BEHIND } from './surfaces.js';

export const CAPACITY = 65536; // segments: a trail of λ = 300 needs about 20,000, or more near the extremal edge
const WIDTH = 3; // CSS pixels
const CASING = 1; // CSS pixels each side
const FADED = 0.15; // the colour left at the tail, as a fraction of the way from BEHIND to the full colour
const FULL = 0.5; // the newest half of the trail is at full colour

// An interleaved, instanced buffer of 6 floats a segment (start, then end), and its two attributes.
function segments(geometry, array, start, end) {
  const buffer = new InstancedInterleavedBuffer(array, 6, 1).setUsage(DynamicDrawUsage);
  geometry.setAttribute(start, new InterleavedBufferAttribute(buffer, 3, 0));
  geometry.setAttribute(end, new InterleavedBufferAttribute(buffer, 3, 3));
  return buffer;
}

// Writes the colour of the angle ψ, faded towards BEHIND by its place along the trail, in linear RGB.
const behind = new Color(BEHIND);
function paint(psi, fraction, out, offset) {
  const k = 3 * wheelIndex(psi);
  const weight = FADED + (1 - FADED) * Math.min(1, fraction / FULL);
  out[offset] = behind.r + weight * (WHEEL[k] - behind.r);
  out[offset + 1] = behind.g + weight * (WHEEL[k + 1] - behind.g);
  out[offset + 2] = behind.b + weight * (WHEEL[k + 2] - behind.b);
}

/** { group, show(path, head, length), hide(), photon }: the trail and the photon, to add to the scene. */
export function createPath() {
  const positions = new Float32Array(6 * CAPACITY);
  const colours = new Float32Array(6 * CAPACITY);
  const geometry = new LineSegmentsGeometry();
  const positionBuffer = segments(geometry, positions, 'instanceStart', 'instanceEnd');
  const colourBuffer = segments(geometry, colours, 'instanceColorStart', 'instanceColorEnd');
  geometry.instanceCount = 0;

  const line = new LineSegments2(geometry, new LineMaterial({
    vertexColors: true, linewidth: WIDTH, fog: true,
    stencilWrite: true, stencilRef: 1, stencilFunc: AlwaysStencilFunc, stencilZPass: ReplaceStencilOp,
  }));
  const casing = new LineSegments2(geometry, new LineMaterial({
    color: 0x07090c, linewidth: WIDTH + 2 * CASING, fog: true,
    stencilWrite: true, stencilRef: 1, stencilFunc: NotEqualStencilFunc,
  }));
  line.renderOrder = 2;
  casing.renderOrder = 3;

  // The photon: a white dot, ringed in black by the inside of a slightly larger sphere.
  const ball = new SphereGeometry(1, 24, 16);
  const ringed = new Mesh(ball, new MeshBasicMaterial({ color: 0x000000, side: BackSide }));
  ringed.scale.setScalar(1.5);
  const photon = new Group().add(new Mesh(ball, new MeshBasicMaterial({ color: 0xffffff })), ringed);

  const group = new Group().add(casing, line, photon);
  for (const mesh of [casing, line]) mesh.frustumCulled = false; // the buffers change every frame
  const head = [0, 0, 0, 0];

  return {
    group,
    photon,
    /** Draws the trail from λ = head − length to the photon at λ = head. */
    show(path, lam, length) {
      const count = fill(path, lam, length, positions, colours, paint, head);
      geometry.instanceCount = count;
      for (const buffer of [positionBuffer, colourBuffer]) {
        buffer.clearUpdateRanges();
        buffer.addUpdateRange(0, 6 * count);
        buffer.needsUpdate = true;
      }
      photon.position.set(head[0], head[1], head[2]);
      group.visible = true;
    },
    hide() {
      group.visible = false;
    },
  };
}
