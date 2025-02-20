// vector ops
const add = (x,y) => x.map((e,i)=>e+y[i])
const sub = (x,y) => x.map((e,i)=>e-y[i])
const mul = (x,s) => x.map(e=>e*s)
const div = (x,s) => x.map(e=>e/s)
const dot = (x,y) => x.map((e,i)=>e*y[i]).reduce((s,t)=>s+t)
const mag = (x) => Math.sqrt(dot(x,x))
const normalize = (x) => div(x,mag(x))
const cross = (x,y) => x.length == 2 ?
  x[0]*y[1]-x[1]*y[0] :
  x.map((e,i)=> x[(i+1)%3]*y[(i+2)%3] - x[(i+2)%3]*y[(i+1)%3])

// matrix ops
const m4row = (m,r) => new m.constructor(4).map((e,i)=>m[r+4*i])
const m4rowdot = (m,r,v) => m[r]*v[0] + m[r+4]*v[1] + m[r+8]*v[2] + m[r+12]*v[3]
const m4col = (m,c) => m.slice(c*4,(c+1)*4)
const m4transpose = (m) => m.map((e,i) => m[((i&3)<<2)+(i>>2)])
const m4mul = (...args) => args.reduce((m1,m2) => {
  if(m2.length == 4) return m2.map((e,i)=>m4rowdot(m1,i,m2)) // m*v
  if(m1.length == 4) return m1.map((e,i)=>m4rowdot(m2,i,m1)) // v*m
  let ans = new m1.constructor(16)
  for(let c=0; c<4; c+=1) for(let r=0; r<4; r+=1)
    ans[r+c*4] = m4rowdot(m1,r,m4col(m2,c))
  return ans // m*m
})

const m3add = (m1, m2) => m1.map((e, i) => e + m2[i]);

const m3identity = () => new Float32Array([
  1, 0, 0,  
  0, 1, 0,  
  0, 0, 1,  
]);

const m3transpose = (m) => {
  return new Float32Array([
    m[0], m[3], m[6], // First row of the transposed matrix
    m[1], m[4], m[7], // Second row of the transposed matrix
    m[2], m[5], m[8]  // Third row of the transposed matrix
  ]);
};

const m4identity = () => new Float32Array([
  1, 0, 0, 0,
  0, 1, 0, 0,
  0, 0, 1, 0,
  0, 0, 0, 1
]);


const m4tom3 = (m4) => {
  return new Float32Array([
    m4[0], m4[1], m4[2], // First row of mat3
    m4[4], m4[5], m4[6], // Second row of mat3
    m4[8], m4[9], m4[10] // Third row of mat3
  ]);
};



// graphics matrices
const m4trans = (dx,dy,dz) => new Float32Array([1,0,0,0, 0,1,0,0, 0,0,1,0, dx,dy,dz,1])
const m4rotX = (ang) => { // around x axis
  let c = Math.cos(ang), s = Math.sin(ang);
  return new Float32Array([1,0,0,0, 0,c,s,0, 0,-s,c,0, 0,0,0,1]);
}
const m4rotY = (ang) => { // around y axis
  let c = Math.cos(ang), s = Math.sin(ang);
  return new Float32Array([c,0,-s,0, 0,1,0,0, s,0,c,0, 0,0,0,1]);
}
const m4rotZ = (ang) => { // around z axis
  let c = Math.cos(ang), s = Math.sin(ang);
  return new Float32Array([c,s,0,0, -s,c,0,0, 0,0,1,0, 0,0,0,1]);
}
const m4fixAxes = (f, up) => { // f to -z, up to near +y
  f = normalize(f)
  let r = normalize(cross(f,up))
  let u = cross(r,f)
  return new Float32Array([
    r[0],u[0],-f[0],0,
    r[1],u[1],-f[1],0,
    r[2],u[2],-f[2],0,
    0,0,0,1
  ])
}
const m4scale = (sx,sy,sz) => new Float32Array([sx,0,0,0, 0,sy,0,0, 0,0,sz,0, 0,0,0,1])
const m4view = (eye, center, up) => m4mul(m4fixAxes(sub(center,eye), up), m4trans(-eye[0],-eye[1],-eye[2]))
const m4perspNegZ = (near, far, fovy, width, height) => {
  let sy = 1/Math.tan(fovy/2);
  let sx = sy*height/width;
  return new Float32Array([sx,0,0,0, 0,sy,0,0, 0,0,-(far+near)/(far-near),-1, 0,0,(2*far*near)/(near-far),0]);
}


const m3mul = (...args) => args.reduce((m1, m2) => {
  if (m2.length === 3) { // m*v (3x3 matrix * 3D vector)
      return m2.map((_, i) => m3rowdot(m1, i, m2));
  }
  if (m1.length === 3) { // v*m (3D vector * 3x3 matrix)
      return m1.map((_, i) => m3rowdot(m2, i, m1));
  }

  // m*m (3x3 matrix * 3x3 matrix)
  let ans = new m1.constructor(9);
  for (let c = 0; c < 3; c++) {
      for (let r = 0; r < 3; r++) {
          ans[r + c * 3] = m3rowdot(m1, r, m3col(m2, c));
      }
  }
  return ans;
});

// Helper function to calculate the dot product of a row and a column
function m3rowdot(matrix, row, vec) {
  return matrix[row] * vec[0] + matrix[row + 3] * vec[1] + matrix[row + 6] * vec[2];
}

// Helper function to extract a column from a 3x3 matrix
function m3col(matrix, col) {
  return [matrix[col * 3], matrix[col * 3 + 1], matrix[col * 3 + 2]];
}


// quaternion
const m4fromQ = (q) => { 
  let n = dot(q,q)
  let s = n ? 2/n : 0
  let xx = s*q[1]*q[1], xy = s*q[1]*q[2], xz = s*q[1]*q[3], xw = s*q[1]*q[0]
  let yy = s*q[2]*q[2], yz = s*q[2]*q[3], yw = s*q[2]*q[0]
  let zz = s*q[3]*q[3], zw = s*q[3]*q[0]
  return new Float32Array([
    1-yy-zz, xy+zw, xz-yw, 0,
    xy-zw, 1-xx-zz, yz+xw, 0,
    xz+yw, yz-xw, 1-xx-yy, 0,
    0,0,0,1,
  ])
}
const m4toQ = (m) => {
  let a00 = m[0], a11 = m[5], a22 = m[10]
  if (a00 + a11 + a22 > 0)
    return normalize([a00+a11+a22+1, m[6]-m[9], m[8]-m[2], m[1]-m[4]])
  if ((a00 > a11) && (a00 > a22))
    return normalize([m[6]-m[9], a00-a11-a22+1, m[1]+m[4], m[8]-m[2]])
  if (a11 > a22)
    return normalize([m[8]-m[2], m[1]+m[4], a11-a00-a22+1, m[6]+m[9]])
  return normalize([m[1]-m[4], m[2]+m[8], m[6]+m[9], a22-a00-a11+1])
}

const m3invert = (m) => {
  // Calculate the determinant of the 3x3 matrix
  const det = m[0] * (m[4] * m[8] - m[5] * m[7]) -
              m[1] * (m[3] * m[8] - m[5] * m[6]) +
              m[2] * (m[3] * m[7] - m[4] * m[6]);

  // If determinant is zero, the matrix is not invertible
  if (det === 0) {
    return null; // No inverse exists
  }

  const invDet = 1 / det;

  // Compute the adjugate matrix (transpose of the cofactor matrix)
  return new Float32Array([
    (m[4] * m[8] - m[5] * m[7]) * invDet,
    (m[2] * m[7] - m[1] * m[8]) * invDet,
    (m[1] * m[5] - m[2] * m[4]) * invDet,

    (m[5] * m[6] - m[3] * m[8]) * invDet,
    (m[0] * m[8] - m[2] * m[6]) * invDet,
    (m[2] * m[3] - m[0] * m[5]) * invDet,

    (m[3] * m[7] - m[4] * m[6]) * invDet,
    (m[1] * m[6] - m[0] * m[7]) * invDet,
    (m[0] * m[4] - m[1] * m[3]) * invDet
  ]);
};


// interpolation
const lerp = (t,p0,p1) => add(mul(p0,1-t), mul(p1,t))
const lbez = (t, ...p) => {
  while(p.length > 1) p = p.slice(1).map((e,i) => lerp(t,p[i],e))
  return p[0]
}
const bezcut = (t, ...p) => {
  let front = [], back = []
  while(p.length > 0) {
    front.push(p[0])
    back.unshift(p[p.length-1])
    p = p.slice(1).map((e,i) => lerp(t,p[i],e))
  }
  return [front, back]
}
const slerp = (t,q0,q1) => {
  let d = dot(q0,q1)
  if (d > 0.9999) return normalize(lerp(t,q0,q1))
  let o = Math.acos(d), den = Math.sin(o)
  return add(mul(q0, Math.sin((1-t)*o)/den), mul(q1, Math.sin(t*o)/den))
}
const qlerp = (t,q0,q1) => {
  let d = dot(q0,q1)
  if (d < 0) { q1 = mul(q1,-1); d = -d; }
  if (d > 0.9999) return normalize(lerp(t,q0,q1))
  let o = Math.acos(d), den = Math.sin(o)
  return add(mul(q0, Math.sin((1-t)*o)/den), mul(q1, Math.sin(t*o)/den))
}
const sbez = (t, ...p) => {
  while(p.length > 1) p = p.slice(1).map((e,i) => slerp(t,p[i],e))
  return p[0]
}
const qbez = (t, ...p) => {
  while(p.length > 1) p = p.slice(1).map((e,i) => qlerp(t,p[i],e))
  return p[0]
}
