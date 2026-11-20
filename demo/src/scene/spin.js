// The hole's spin a: an arrow below the orbit's spheroid, around a short axle on the z axis, the way the hole turns:
// towards −φ (docs/equations.md, "Signs"), which is clockwise seen from above. It is below because most of the
// presets' paths are on the upper half. It hangs far enough down to clear the spheroid's outline from every
// preset's viewpoint, which look down at up to 36°; from higher up, the spheroid hides it. The spin b turns ψ,
// which the 3D view shows only by colour, so its arrow goes round the colour wheel instead.
//
// The arrow is a muted grey-blue with a dark outline, so it reads on the dark spheroid and on either theme's page
// without drawing the eye from the path. The outline is a slightly fatter copy of each part with only its inside
// faces drawn. The arrow turns with the camera about the axle, so its head is always at the front left: the view
// turns, and the arrow stays where it is on screen.

import { BackSide, ConeGeometry, CylinderGeometry, Group, Mesh, MeshBasicMaterial, SphereGeometry, TorusGeometry } from 'three';

import { SURFACE } from './surfaces.js';

/** The spin arrows' colours, here and round the colour wheel. */
export const ARROW = { core: SURFACE.rim, edge: SURFACE.backdrop };

const RADIUS = 0.26; // the arrow's, about the axle
const DEPTH = 0.6; // the arrow's, below the bottom of the spheroid: the Equal spins preset needs 0.6, the rest less
const AXLE = [0.47, 0.73]; // where the axle starts and stops, below the bottom of the spheroid
const TUBE = 0.013; // the radius of the arrow's and the axle's stroke
const EDGE = 0.01; // the outline's width
const ARC = (290 * Math.PI) / 180; // the arrow's length, in angle
const HEAD = [0.047, 0.11]; // the head's radius and length
const FRONT = (-50 * Math.PI) / 180; // where the head is, from the side facing the camera: front left

const core = new MeshBasicMaterial({ color: ARROW.core, fog: false });
const edge = new MeshBasicMaterial({ color: ARROW.edge, side: BackSide, fog: false });

// A part and its outline, made by `shape(grow)` with its sizes grown by `grow`.
function outlined(shape) {
  return new Group().add(new Mesh(shape(0), core), new Mesh(shape(EDGE), edge));
}

/**
 * The axle and arrow: { group, setShape(up, spinning), face(camera) }. setShape puts them below a spheroid whose
 * semi-axis along z is `up`, and shows them only while the hole spins in the a-plane. face(camera) turns the
 * arrow's head to the front.
 */
export function createSpin() {
  const length = AXLE[1] - AXLE[0];
  const axle = outlined((grow) => new CylinderGeometry(TUBE + grow, TUBE + grow, length + 2 * grow, 12).rotateX(Math.PI / 2));
  axle.position.z = -AXLE[0] - length / 2;

  // The torus starts at angle 0 and runs counter-clockwise, so the head goes at its start, pointing clockwise (−y).
  const arrow = new Group().add(
    outlined((grow) => new TorusGeometry(RADIUS, TUBE + grow, 12, 96, ARC)),
    outlined((grow) => new SphereGeometry(TUBE + grow, 12, 8).translate(RADIUS * Math.cos(ARC), RADIUS * Math.sin(ARC), 0)),
    outlined((grow) => new ConeGeometry(HEAD[0] + 1.5 * grow, HEAD[1] + 3 * grow, 24).rotateZ(Math.PI).translate(RADIUS, -HEAD[1] / 2, 0)),
  );
  arrow.position.z = -DEPTH;

  const group = new Group().add(axle, arrow);
  return {
    group,
    setShape(up, spinning) {
      group.visible = spinning;
      group.position.z = -up;
    },
    face(camera) {
      arrow.rotation.z = Math.atan2(camera.position.y, camera.position.x) + FRONT;
    },
  };
}
