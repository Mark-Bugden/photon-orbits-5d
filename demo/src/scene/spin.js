// The hole's spin a: an arrow above the orbit's spheroid, around a short axle, the way the hole turns: towards −φ
// (docs/equations.md, "Signs"), which is clockwise seen from above. The spin b turns ψ, which the 3D view shows
// only by colour. The axle stops short of the spheroid, so it never hides a path over the pole.
//
// The arrow is pale with a dark outline, so it reads on the dark spheroid and on either theme's page. The outline
// is a slightly fatter copy of each part with only its inside faces drawn. The arrow turns with the camera about
// the axle, so its head is always at the front left: the view turns, and the arrow stays where it is on screen.

import { BackSide, ConeGeometry, CylinderGeometry, Group, Mesh, MeshBasicMaterial, SphereGeometry, TorusGeometry } from 'three';

import { SURFACE } from './surfaces.js';

const RADIUS = 0.3; // the arrow's, about the axle
const HEIGHT = 0.25; // the arrow's, above the top of the spheroid
const AXLE = [0.12, 0.42]; // where the axle starts and stops, above the top of the spheroid
const TUBE = 0.016; // the radius of the arrow's and the axle's stroke
const EDGE = 0.012; // the outline's width
const ARC = (290 * Math.PI) / 180; // the arrow's length, in angle
const HEAD = [0.055, 0.13]; // the head's radius and length
const FRONT = (-50 * Math.PI) / 180; // where the head is, from the side facing the camera: front left

const pale = new MeshBasicMaterial({ color: SURFACE.ring, fog: false });
const dark = new MeshBasicMaterial({ color: SURFACE.backdrop, side: BackSide, fog: false });

// A part and its outline, made by `shape(grow)` with its sizes grown by `grow`.
function outlined(shape) {
  return new Group().add(new Mesh(shape(0), pale), new Mesh(shape(EDGE), dark));
}

/**
 * The axle and arrow: { group, setShape(top, spinning), face(camera) }. setShape puts them on top of a spheroid
 * whose top is at z = top, and shows them only while the hole spins in the a-plane. face(camera) turns the arrow's
 * head to the front.
 */
export function createSpin() {
  const length = AXLE[1] - AXLE[0];
  const axle = outlined((grow) => new CylinderGeometry(TUBE + grow, TUBE + grow, length + 2 * grow, 12).rotateX(Math.PI / 2));
  axle.position.z = AXLE[0] + length / 2;

  // The torus starts at angle 0 and runs counter-clockwise, so the head goes at its start, pointing clockwise (−y).
  const arrow = new Group().add(
    outlined((grow) => new TorusGeometry(RADIUS, TUBE + grow, 12, 96, ARC)),
    outlined((grow) => new SphereGeometry(TUBE + grow, 12, 8).translate(RADIUS * Math.cos(ARC), RADIUS * Math.sin(ARC), 0)),
    outlined((grow) => new ConeGeometry(HEAD[0] + 1.5 * grow, HEAD[1] + 3 * grow, 24).rotateZ(Math.PI).translate(RADIUS, -HEAD[1] / 2, 0)),
  );
  arrow.position.z = HEIGHT;

  const group = new Group().add(axle, arrow);
  return {
    group,
    setShape(top, spinning) {
      group.visible = spinning;
      group.position.z = top;
    },
    face(camera) {
      arrow.rotation.z = Math.atan2(camera.position.y, camera.position.x) + FRONT;
    },
  };
}
