// The 3D view: the spheroids, the trail and the photon, a camera that turns about the hole, and its buttons.
// Drag to turn it, or use the arrow keys; zoom with the buttons, by pinching, or with + and −. The mouse wheel
// scrolls the page as usual.

import { Color, Fog, PerspectiveCamera, Scene, Vector3, WebGLRenderer } from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

import { createPath } from './path.js';
import { BEHIND, createSurfaces } from './surfaces.js';

const FOV = 30; // degrees, vertical
const VIEWPOINT = new Vector3(1.3, -2.4, 2).normalize(); // Mathematica's default, until a preset sets another
const MARGIN = 1.12; // room around the spheroid when the view fits it
const ZOOM = 0.8; // each press of a zoom button
const TURN = Math.PI / 24; // each press of an arrow key: 7.5°

const ICONS = {
  in: '<path d="M12 5v14M5 12h14"/>',
  out: '<path d="M5 12h14"/>',
  fit: '<path d="M4 9V4h5M20 9V4h-5M4 15v5h5M20 15v5h-5"/><circle cx="12" cy="12" r="2.5"/>',
};

function button(label, icon, action) {
  const element = document.createElement('button');
  element.type = 'button';
  element.className = 'dk-button po-view-button';
  element.title = label;
  element.setAttribute('aria-label', label);
  element.innerHTML = `<svg viewBox="0 0 24 24" width="18" height="18" aria-hidden="true" focusable="false" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">${ICONS[icon]}</svg>`;
  element.addEventListener('click', action);
  return element;
}

/**
 * Creates the view in `stage`. Throws if WebGL can't start. onChange() is called whenever the camera moves, so the
 * caller can draw again. Returns { setSize, setShape, fit, setBackground, showTrail, hideTrail, render }.
 */
export function createView(stage, { onChange, reducedMotion = false }) {
  // The canvas is cleared to the stage's colour (setBackground) rather than left clear, so the page can't show
  // through anything drawn with partial transparency.
  const renderer = new WebGLRenderer({ antialias: true });
  const canvas = renderer.domElement;
  canvas.className = 'po-canvas';
  canvas.tabIndex = 0;
  canvas.setAttribute('aria-label', '3D view of the orbit. Drag, or use the arrow keys, to turn it; + and − zoom.');

  const scene = new Scene();
  scene.fog = new Fog(BEHIND, 1, 10);
  const surfaces = createSurfaces();
  const path = createPath();
  scene.add(surfaces.horizon, surfaces.orbit, path.group);

  const camera = new PerspectiveCamera(FOV, 1, 0.01, 200);
  camera.up.set(0, 0, 1); // before the controls, which turn about camera.up
  const controls = new OrbitControls(camera, canvas);
  controls.enablePan = false;
  controls.enableDamping = !reducedMotion;
  controls.dampingFactor = 0.12;
  controls.rotateSpeed = 0.8;
  controls.minDistance = 1;
  controls.maxDistance = 60;
  controls.addEventListener('change', onChange);

  // The wheel scrolls the page: stop it before it reaches the controls, without preventing the scroll.
  stage.addEventListener('wheel', (event) => event.stopPropagation(), { capture: true });

  const zoomIn = () => controls.dollyIn(ZOOM);
  const zoomOut = () => controls.dollyOut(ZOOM);
  const keys = {
    ArrowLeft: () => controls.rotateLeft(TURN),
    ArrowRight: () => controls.rotateLeft(-TURN),
    ArrowUp: () => controls.rotateUp(TURN),
    ArrowDown: () => controls.rotateUp(-TURN),
    '+': zoomIn,
    '=': zoomIn,
    '-': zoomOut,
    _: zoomOut,
  };
  canvas.addEventListener('keydown', (event) => {
    const action = keys[event.key];
    if (!action || event.altKey || event.ctrlKey || event.metaKey) return;
    event.preventDefault();
    action();
  });

  let radius = 1; // of the orbit's spheroid, or the horizon's if there is no orbit: what the view fits
  let size = { width: 0, height: 0, pixelRatio: 0 };
  const home = VIEWPOINT.clone(); // the direction "Reset view" looks from

  // Turns the camera to look from home, or from [azimuth, elevation] in degrees, which becomes home, at the
  // distance that fits the spheroid.
  function fit(view) {
    if (Array.isArray(view)) {
      const [azimuth, elevation] = view.map((degrees) => (degrees * Math.PI) / 180);
      home.set(Math.cos(elevation) * Math.cos(azimuth), Math.cos(elevation) * Math.sin(azimuth), Math.sin(elevation));
    }
    const half = (FOV * Math.PI) / 360;
    const across = Math.atan(Math.tan(half) * camera.aspect);
    const damping = controls.enableDamping;
    controls.enableDamping = false;
    controls.update(); // finishes any turn still easing out, so it can't carry on after the camera moves
    camera.position.copy(home).multiplyScalar((MARGIN * radius) / Math.sin(Math.min(half, across)));
    controls.target.set(0, 0, 0);
    controls.update();
    controls.enableDamping = damping;
    onChange();
  }

  const buttons = document.createElement('div');
  buttons.className = 'po-view-buttons';
  buttons.append(button('Zoom in', 'in', zoomIn), button('Zoom out', 'out', zoomOut), button('Reset view', 'fit', fit));
  stage.append(canvas, buttons);

  return {
    /** Matches the drawing buffer to the stage. Cheap when nothing changed. */
    setSize(width, height, pixelRatio) {
      if (width === size.width && height === size.height && pixelRatio === size.pixelRatio) return;
      size = { width, height, pixelRatio };
      renderer.setPixelRatio(Math.min(pixelRatio, 2));
      renderer.setSize(Math.max(1, width), Math.max(1, height), false);
      camera.aspect = Math.max(1, width) / Math.max(1, height);
      camera.updateProjectionMatrix();
    },
    /** Shapes the spheroids, each given as its semi-axes [across, up], or null to hide the orbit's. */
    setShape(orbitAxes, horizonAxes) {
      surfaces.setShape(orbitAxes, horizonAxes);
      if (orbitAxes) path.setSpheroid(orbitAxes);
      radius = Math.max(...(orbitAxes ?? horizonAxes));
      path.photon.scale.setScalar(0.022 * radius);
    },
    /**
     * Turns the camera to look from [azimuth, elevation], in degrees with z up, at the distance that fits the
     * current spheroid; "Reset view" then comes back here. Without a view, it looks from the last one given.
     */
    fit,
    /** Paints the background around the spheroid, e.g. the stage's colour for the theme. */
    setBackground(colour) {
      renderer.setClearColor(new Color(colour));
    },
    showTrail: path.show,
    hideTrail: path.hide,
    render() {
      controls.update(); // settles the damping; calls onChange while the camera still moves
      // Faint fog: none on the near side of the spheroid, about a third of the way to BEHIND on the far side.
      const distance = camera.position.distanceTo(controls.target);
      scene.fog.near = distance - 0.2 * radius;
      scene.fog.far = distance + 3.5 * radius;
      // Near and far planes close around the scene, for the depth buffer's precision.
      camera.near = Math.max(0.01, distance - 3 * radius);
      camera.far = distance + 3 * radius;
      camera.updateProjectionMatrix();
      renderer.render(scene, camera);
    },
  };
}
