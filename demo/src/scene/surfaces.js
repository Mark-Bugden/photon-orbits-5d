// The orbit's spheroid and its equator ring, and the horizon inside, after the paper's figures. They are dark in
// both themes, so the path's colours read the same in both, and nothing behind the spheroid shows through it.
//
// The path lies on the spheroid, so the spheroid is drawn as two shells either side of it:
// - a backdrop just outside, inside faces only and opaque: the far half of the path lies on it;
// - a translucent shell just inside, outside faces only and drawn last: the near half of the path lies over it,
//   and it dims the far half, which shows through it.

import { BackSide, Color, Group, Mesh, MeshBasicMaterial, ShaderMaterial, SphereGeometry } from 'three';
import { LineMaterial } from 'three/addons/lines/LineMaterial.js';
import { LineSegments2 } from 'three/addons/lines/LineSegments2.js';
import { LineSegmentsGeometry } from 'three/addons/lines/LineSegmentsGeometry.js';

export const SURFACE = {
  backdrop: '#11151d', // the inside of the far half
  shell: '#3a4558', // the near half, shaded
  rim: '#7b8aa5', // brightens the spheroid's edge, so its outline reads on the dark page
  alpha: 0.62, // the near half's opacity
  horizon: '#000000',
  horizonRim: '#2a3240',
  ring: '#dde2ea', // the equator
};

// What the trail fades towards and the fog fades to: about the spheroid's colour where the path lies on it.
export const BEHIND = '#1e2531';

const OUTSIDE = 1.006; // the backdrop's scale, and the shell's, either side of the path
const INSIDE = 0.994;

const vertexShader = /* glsl */ `
  varying vec3 vNormal;
  varying vec3 vView;
  void main() {
    vec4 position4 = modelViewMatrix * vec4(position, 1.0);
    vNormal = normalize(normalMatrix * normal);
    vView = -position4.xyz;
    gl_Position = projectionMatrix * position4;
  }
`;

// Matte shading from a light above and to the left of the camera, and a rim that brightens towards the outline.
const fragmentShader = /* glsl */ `
  uniform vec3 base;
  uniform vec3 rim;
  uniform float alpha;
  varying vec3 vNormal;
  varying vec3 vView;
  void main() {
    vec3 normal = normalize(vNormal);
    float facing = clamp(dot(normal, normalize(vView)), 0.0, 1.0);
    float light = 0.55 + 0.45 * max(dot(normal, normalize(vec3(-0.35, 0.55, 0.75))), 0.0);
    gl_FragColor = vec4(base * light + rim * pow(1.0 - facing, 3.0), alpha);
    #include <colorspace_fragment>
  }
`;

function shaded(base, rim, alpha) {
  return new ShaderMaterial({
    uniforms: { base: { value: new Color(base) }, rim: { value: new Color(rim) }, alpha: { value: alpha } },
    vertexShader,
    fragmentShader,
    transparent: alpha < 1,
    depthWrite: alpha >= 1,
  });
}

// A unit sphere with its poles on the z axis, which is up in the scene.
const sphere = new SphereGeometry(1, 96, 64).rotateX(Math.PI / 2);

function equator(segments = 256) {
  const points = [];
  for (let k = 0; k < segments; k++) {
    const [from, to] = [k, k + 1].map((j) => (2 * Math.PI * j) / segments);
    points.push(Math.cos(from), Math.sin(from), 0, Math.cos(to), Math.sin(to), 0);
  }
  return new LineSegmentsGeometry().setPositions(points);
}

/**
 * The spheroids: { orbit, horizon }, two groups to add to the scene, and setShape(orbit, horizon), which scales
 * each to its semi-axes [across, up], or hides it for null.
 */
export function createSurfaces() {
  const backdrop = new Mesh(sphere, new MeshBasicMaterial({ color: SURFACE.backdrop, side: BackSide }));
  backdrop.scale.setScalar(OUTSIDE);
  const shell = new Mesh(sphere, shaded(SURFACE.shell, SURFACE.rim, SURFACE.alpha));
  shell.scale.setScalar(INSIDE);
  shell.renderOrder = 10;
  const ring = new LineSegments2(equator(), new LineMaterial({ color: SURFACE.ring, linewidth: 1.5, fog: true, alphaToCoverage: true }));
  ring.scale.setScalar(1.002);
  const orbit = new Group().add(backdrop, ring, shell);
  const horizon = new Group().add(new Mesh(sphere, shaded(SURFACE.horizon, SURFACE.horizonRim, 1)));

  const place = (group, axes) => {
    group.visible = axes !== null;
    if (axes) group.scale.set(axes[0], axes[0], axes[1]);
  };
  return {
    orbit,
    horizon,
    setShape(orbitAxes, horizonAxes) {
      place(orbit, orbitAxes);
      place(horizon, horizonAxes);
    },
  };
}
