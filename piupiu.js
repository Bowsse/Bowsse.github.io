! function ()
{
	"use strict";

	var screen = ge1doot.screen.init("screen", function ()
	{
		PHY2D.zoom = 1;
		//PHY2D.deleteStatic();
		PHY2D.rectangle(screen.width / 2, screen.height + 50, screen.width, 100, 0, 0);
		PHY2D.rectangle(screen.width / 2, -screen.height / 2, screen.width, 100, 0, 0);
		PHY2D.rectangle(-50, screen.height / 2, 100, screen.height * 3, 0, 0);
		PHY2D.rectangle(screen.width + 50, screen.height / 2 , 100, screen.height * 3, 0, 0);
		PHY2D.zoom = screen.width / 1920;
	}, false);
	var ctx = screen.ctx;
	if (!ctx.setLineDash)
	{
		ctx.setLineDash = function () {}
	}
	var pointer = screen.pointer.init(
	{
		up: function ()
		{
			PHY2D.stopManipulate();
		}
	});

	// Main PHY2D code
	
	var PHY2D = {};
	
	(function(ctx, pointer, kGravity, kTimeStep, kFriction, kMotionAABB, kImgPath) {

		var objects = [];

		function Vector(x, y)
		{
			this.x = x || 0.0;
			this.y = y || 0.0;
		}
		Vector.prototype = {
			set: function (x, y)
			{
				this.x = x;
				this.y = y;
				return this;
			},
			dot: function (v)
			{
				return this.x * v.x + this.y * v.y;
			},
			lenSqr: function ()
			{
				return this.x * this.x + this.y * this.y;
			},
			len: function ()
			{
				return Math.sqrt(this.x * this.x + this.y * this.y);
			},
			transform: function (v, m)
			{
				this.x = m.cos * v.x - m.sin * v.y + m.pos.x;
				this.y = m.sin * v.x + m.cos * v.y + m.pos.y;
				return this;
			},
			rotate: function (v, m)
			{
				this.x = m.cos * v.x - m.sin * v.y;
				this.y = m.sin * v.x + m.cos * v.y;
				return this;
			},
			normal: function (a, b)
			{
				var x = a.x - b.x,
					y = a.y - b.y,
					len = Math.sqrt(x * x + y * y);
				this.x = -y / len;
				this.y = x / len;
				return this;
			},
			project: function (a, b, n)
			{
				var x = a.x - b.x,
					y = a.y - b.y,
					len = Math.sqrt(x * x + y * y);
				return (-y / len) * n.x + (x / len) * n.y;
			},
			addScale: function (v1, v2, s)
			{
				this.x = v1.x + (v2.x * s);
				this.y = v1.y + (v2.y * s);
				return this;
			},
			subScale: function (v1, v2, s)
			{
				this.x = v1.x - (v2.x * s);
				this.y = v1.y - (v2.y * s);
				return this;
			},
			add: function (v1, v2)
			{
				this.x = v1.x + v2.x;
				this.y = v1.y + v2.y;
				return this;
			},
			sub: function (v1, v2)
			{
				this.x = v1.x - v2.x;
				this.y = v1.y - v2.y;
				return this;
			},
			scale: function (v1, s)
			{
				this.x = v1.x * s;
				this.y = v1.y * s;
				return this;
			},
			perp: function ()
			{
				var x = this.x;
				this.x = -this.y;
				this.y = x;
				return this;
			},
			inv: function (v1)
			{
				this.x = -v1.x;
				this.y = -v1.y;
				return this;
			},
			clamp: function (v, min, max)
			{
				if (v > max) v = max;
				else if (v < min) v = min;
				return v;
			},
			rotateIntoSpaceOf: function (a, m)
			{
				var dx = -a.x,
					dy = -a.y;
				this.x = dx * m.cos + dy * m.sin;
				this.y = dx * -m.sin + dy * m.cos;
				return this;
			},
			// SIMD array vectors
			array: function (n, values)
			{
				var array = new Array(n);
				array.min = new Vector();
				array.max = new Vector();
				for (var i = 0; i < n; i++)
				{
					array[i] = new Vector(
						values ? values[i * 2 + 0] : 0.0,
						values ? values[i * 2 + 1] : 0.0
					);
				}
				array.transform = function (v, m)
				{
					for (var i = 0, len = this.length; i < len; i++)
					{
						var vi = v[i],
							elem = this[i];
						var x = m.cos * vi.x - m.sin * vi.y + m.pos.x;
						var y = m.sin * vi.x + m.cos * vi.y + m.pos.y;
						if (x < this.min.x) this.min.x = x;
						if (y < this.min.y) this.min.y = y;
						if (x > this.max.x) this.max.x = x;
						if (y > this.max.y) this.max.y = y;
						elem.x = x;
						elem.y = y;
					}
					return this;
				}
				array.rotate = function (v, m)
				{
					for (var i = 0, len = this.length; i < len; i++)
					{
						var vi = v[i],
							elem = this[i];
						elem.x = m.cos * vi.x - m.sin * vi.y;
						elem.y = m.sin * vi.x + m.cos * vi.y;
					}
					return this;
				}
				array.resetMinmax = function ()
				{
					this.min.x = 10000000000.0;
					this.min.y = 10000000000.0;
					this.max.x = -10000000000.0;
					this.max.y = -10000000000.0;
				}
				array.normal = function (points)
				{
					for (var i = 0; i < this.length; i++)
					{
						this[i].normal(
							points[(i + 1) % this.length],
							points[i]
						);
					}
					return this;
				}
				return array;
			}
		}
		// Matrix container (but it's not a matrix)
		function Matrix()
		{
			this.cos = 0.0;
			this.sin = 0.0;
			this.pos = new Vector();
			this.ang = 0.0;
		}
		Matrix.prototype = {
			set: function (a, x, y)
			{
				this.cos = Math.cos(a);
				this.sin = Math.sin(a);
				this.ang = a;
				this.pos.x = x;
				this.pos.y = y;
				return this;
			},
			copy: function (matrix)
			{
				this.cos = matrix.cos;
				this.sin = matrix.sin;
				this.ang = matrix.ang;
				this.pos.x = matrix.pos.x;
				this.pos.y = matrix.pos.y;
				return this;
			},
			integrate: function (va, vx, vy, kTimeStep)
			{
				this.pos.x += vx * kTimeStep;
				this.pos.y += vy * kTimeStep;
				this.ang += va * kTimeStep;
				this.cos = Math.cos(this.ang);
				this.sin = Math.sin(this.ang);
				return this;
			}
		}
		
		var v0 = new Vector();
		var v1 = new Vector();
		var v2 = new Vector();
		var v3 = new Vector();
		var v4 = new Vector();
		var v5 = new Vector();

		// contacts list
		var manipulate = null;
		var joints = [];
		var contacts = [];
		contacts.index = 0;
		contacts.create = function (A, B, pa, pb, nx, ny)
		{
			if (!this[this.index]) this[this.index] = new Contact();
			this[this.index++].set(A, B, pa, pb, nx, ny);
		}

		// AABB container constructor
		function AABB()
		{
			this.x = 0.0;
			this.y = 0.0;
			this.w = 0.0;
			this.h = 0.0;
		}
		

		// Polygon constructor
		function Polygon(x, y, w, h, vertices, mass, angle, img)
		{
			if (img)
			{
				var image = new Image();
				image.src = kImgPath + img.src;
				this.img = {
					w: (img.w ? img.w*PHY2D.zoom : w),
					h: (img.h ? img.h*PHY2D.zoom : h),
					x: (img.x ? img.x*PHY2D.zoom : 0),
					y: (img.y ? img.y*PHY2D.zoom : 0),
					elem: image
				}
			}
			else this.img = null;
			this.vel = new Vector();
			this.angularVel = 0.0;
			this.invMass = mass ? 1 / mass : 0;
			this.matrix = new Matrix().set(angle, x, y);
			kMotionAABB && (this.matrixNextFrame = new Matrix());
			this.aabb = new AABB();
			this.static = false;
			this.length = (vertices.length / 2) | 0;

			// vertices
			this.localSpacePoints = new Vector().array(this.length, vertices);
			this.localSpaceNormals = new Vector().array(this.length).normal(this.localSpacePoints);
			this.worldSpaceNormals = new Vector().array(this.length);
			this.worldSpacePoints = new Vector().array(this.length);

			// calculate inverse inertia tensor
			this.invI = (this.invMass > 0) ? 1 / ((1 / this.invMass) * (w * w + h * h) / 3) : 0

			// contact points
			this.c1 = new Vector();
			this.c0 = new Vector();

			// add rigid body
			objects.push(this);
			
		}

		Polygon.prototype = {

			// aabb motion box

			motionAABB: function ()
			{
				this.worldSpacePoints.resetMinmax();
				kMotionAABB && this.worldSpacePoints.transform(this.localSpacePoints, this.matrixNextFrame);
				this.worldSpacePoints.transform(this.localSpacePoints, this.matrix);
				this.worldSpaceNormals.rotate(this.localSpaceNormals, this.matrix);
				var min = this.worldSpacePoints.min;
				var max = this.worldSpacePoints.max;
				this.aabb.x = (min.x + max.x) * 0.5;
				this.aabb.y = (min.y + max.y) * 0.5;
				this.aabb.w = (max.x - min.x) * 0.5;
				this.aabb.h = (max.y - min.y) * 0.5;
			},
			
			isPointInPoly: function (x, y)
			{
		
			var c = false;
				for (var p = this.worldSpacePoints, i = -1, l = this.length, j = l - 1; ++i < l; j = i)
				{
					((p[i].y <= y && y < p[j].y) || (p[j].y <= y && y < p[i].y))
					&& (x <= (p[j].x - p[i].x) * (y - p[i].y) / (p[j].y - p[i].y) + p[i].x)
					&& (c = !c);
				}
				return c;
			},

			// contact points

			contact: function (that)
			{
				var face, vertex, vertexRect, faceRect, fp, va, vb, vc, nx, ny, wsN, wdV0, wdV1, wsV0, wsV1;
				
				// generate contacts for this pair
				mostSeparated.set(100000, -1, -1, 0, 100000);
				mostPenetrating.set(-100000, -1, -1, 0, 100000);

				// face of A, vertices of B
				this.featurePairJudgement(that, 2);

				// faces of B, vertices of A
				that.featurePairJudgement(this, 1);

				if (mostSeparated.dist > 0 && mostSeparated.fpc !== 0)
				{
					// objects are separated
					face = mostSeparated.edge;
					vertex = mostSeparated.closestI;
					fp = mostSeparated.fpc;
				}
				else if (mostPenetrating.dist <= 0)
				{
					// objects are penetrating
					face = mostPenetrating.edge;
					vertex = mostPenetrating.closestI;
					fp = mostPenetrating.fpc;
				}

				if (fp === 1) vertexRect = this, faceRect = that;
				else vertexRect = that, faceRect = this;

				// world space vertex
				wsN = faceRect.worldSpaceNormals[face];

				// other vertex adjacent which makes most parallel normal with the collision normal
				va = vertexRect.worldSpacePoints[(vertex - 1 + vertexRect.length) % vertexRect.length];
				vb = vertexRect.worldSpacePoints[vertex];
				vc = vertexRect.worldSpacePoints[(vertex + 1) % vertexRect.length];

				if (v0.project(vb, va, wsN) < v1.project(vc, vb, wsN))
				{
					wdV0 = va;
					wdV1 = vb;
				}
				else
				{
					wdV0 = vb;
					wdV1 = vc;
				}

				// world space edge
				wsV0 = faceRect.worldSpacePoints[face];
				wsV1 = faceRect.worldSpacePoints[(face + 1) % faceRect.length];

				// form contact
				if (fp === 1)
				{
					// project vertex onto edge
					this.projectPointOntoEdge(wsV0, wsV1, wdV0, wdV1);
					that.projectPointOntoEdge(wdV1, wdV0, wsV0, wsV1);
					// normal is negated because it always needs to point from A->B
					nx = -wsN.x;
					ny = -wsN.y;
				}
				else
				{
					this.projectPointOntoEdge(wdV1, wdV0, wsV0, wsV1);
					that.projectPointOntoEdge(wsV0, wsV1, wdV0, wdV1);
					nx = wsN.x;
					ny = wsN.y;
				}

				// create contacts
				contacts.create(this, that, this.c0, that.c0, nx, ny);
				contacts.create(this, that, this.c1, that.c1, nx, ny);
			},

			featurePairJudgement: function (that, fpc)
			{
				var wsN, closestI, closest, dist;
				for (var edge = 0; edge < this.length; edge++)
				{
					// get support vertices
					wsN = this.worldSpaceNormals[edge];

					// rotate into RigidBody space
					v5.rotateIntoSpaceOf(wsN, that.matrix);

					var closestI = -1,
						closestD = -100000;

					// Get the vertex most in the direction of the given vector
					for (var i = 0; i < that.length; i++)
					{
						var d = v5.dot(that.localSpacePoints[i]);
						if (d > closestD)
						{
							closestD = d;
							closestI = i;
						}
					}

					var closest = that.worldSpacePoints[closestI];
					v0.sub(closest, this.worldSpacePoints[edge]);

					// distance from origin to face	
					var dist = v0.dot(wsN);

					if (dist > 0)
					{
						// recompute distance to clamped edge
						v1.sub(closest, this.worldSpacePoints[(edge + 1) % this.length]);

						// project onto minkowski edge
						var l = v2.sub(v1, v0).lenSqr() + 0.0000001;
						dist = this.c0.addScale(v0, v2, v3.clamp(v3.inv(v0).dot(v2) / l, 0, 1)).lenSqr();

						// track separation
						if (dist < mostSeparated.dist)
						{
							mostSeparated.set(dist, closestI, edge, fpc);
						}
					}
					else
					{
						// track penetration
						if (dist > mostPenetrating.dist)
						{
							mostPenetrating.set(dist, closestI, edge, fpc);
						}
					}
				}
			},

			projectPointOntoEdge: function (p0, p1, e0, e1)
			{
				var l = v2.sub(e1, e0).lenSqr() + 0.0000001;
				this.c0.addScale(e0, v2, v3.clamp(v3.sub(p0, e0).dot(v2) / l, 0, 1));
				this.c1.addScale(e0, v2, v3.clamp(v3.sub(p1, e0).dot(v2) / l, 0, 1));
			},

			// integration

			integrate: function ()
			{
				// gravity
				if (this.invMass > 0) this.vel.y += kGravity;

				// update position
				this.matrix.integrate(this.angularVel, this.vel.x, this.vel.y, kTimeStep);
				kMotionAABB && this.matrixNextFrame.copy(this.matrix).integrate(this.angularVel, this.vel.x, this.vel.y, kTimeStep);

				// compute motion AABB
				if (!this.static) this.motionAABB();
				else
				{
					if (this.invMass === 0)
					{
						this.static = true;
						this.motionAABB();
					}
				}
				
				// manipulate
				if (pointer.active && !manipulate && this.invMass)
				{
					if (this.isPointInPoly(pointer.pos.x, pointer.pos.y))
					{
						manipulate = new Manipulate(this);
						manipulate.solve();
					}
				}
			},

			draw: function ()
			{
				if (this.img)
				{
					// draw image
					var m = this.matrix, img = this.img;
					ctx.save();
					ctx.translate(m.pos.x, m.pos.y);
					ctx.rotate(m.ang);
					ctx.drawImage(img.elem, -img.w * 0.5 + img.x, -img.h * 0.5 + img.y, img.w, img.h);
					ctx.restore();
					

				}
			}
		}
		
		// feature pair container

		function FeaturePair()
		{
			this.dist = 0;
			this.closestI = 0;
			this.edge = 0;
			this.fpc = 0;
		}
		FeaturePair.prototype.set = function (dist, closestI, edge, fpc)
		{
			this.dist = dist;
			this.closestI = closestI;
			this.edge = edge;
			this.fpc = fpc;
		}
		var mostSeparated = new FeaturePair();
		var mostPenetrating = new FeaturePair();

		// contacts constructor

		function Contact()
		{
			this.a           = null;
			this.b           = null;
			this.normal      = new Vector();
			this.normalPerp  = new Vector();
			this.ra          = new Vector();
			this.rb          = new Vector();
			this.impulse     = new Vector();
			this.dv          = new Vector();
			this.dist        = 0;
			this.impulseN    = 0;
			this.impulseT    = 0;
			this.invDenom    = 0;
			this.invDenomTan = 0;
		}
		Contact.prototype = {

			// reusing existing contact objects

			set: function (A, B, pa, pb, nx, ny)
			{

				var ran, rbn;
				this.a = A;
				this.b = B;
				this.normal.set(nx, ny);
				this.normalPerp.set(-ny, nx);
				this.dist = v1.sub(pb, pa).dot(this.normal);
				this.impulseN = 0;
				this.impulseT = 0;

				// calculate radius arms
				this.ra.sub(pa, A.matrix.pos).perp();
				this.rb.sub(pb, B.matrix.pos).perp();

				// compute denominator in impulse equation
				ran = this.ra.dot(this.normal);
				rbn = this.rb.dot(this.normal);
				this.invDenom = 1 / (A.invMass + B.invMass + (ran * ran * A.invI) + (rbn * rbn * B.invI));
				ran = this.ra.dot(this.normalPerp);
				rbn = this.rb.dot(this.normalPerp);
				this.invDenomTan = 1 / (A.invMass + B.invMass + (ran * ran * A.invI) + (rbn * rbn * B.invI));
			},

			applyImpulse: function ()
			{
				// linear
				this.a.vel.addScale(this.a.vel, this.impulse, this.a.invMass);
				this.b.vel.subScale(this.b.vel, this.impulse, this.b.invMass);
				// angular
				this.a.angularVel += this.impulse.dot(this.ra) * this.a.invI;
				this.b.angularVel -= this.impulse.dot(this.rb) * this.b.invI;
			},

			// speculative contact solver

			solve: function ()
			{
				var newImpulse, absMag;

				// get all of relative normal velocity
				this.dv.sub(
					v1.addScale(this.b.vel, this.rb, this.b.angularVel),
					v2.addScale(this.a.vel, this.ra, this.a.angularVel)
				);

				// accumulated impulse
				newImpulse = (this.dv.dot(this.normal) + this.dist / kTimeStep) * this.invDenom + this.impulseN;

				// push only
				if (newImpulse > 0) newImpulse = 0;

				// apply impulse
				this.impulse.scale(this.normal, newImpulse - this.impulseN)
				this.applyImpulse();
				this.impulseN = newImpulse;

				// friction
				absMag = Math.abs(this.impulseN) * kFriction;
				newImpulse = this.impulse.clamp(this.dv.dot(this.normalPerp) * this.invDenomTan + this.impulseT, -absMag, absMag);

				// apply friction impulse
				this.impulse.scale(this.normalPerp, newImpulse - this.impulseT)
				this.applyImpulse();
				this.impulseT = newImpulse;
			}
		}
		
		// join constructor
		
		function Constraint (A, va, B, vb)
		{
			this.a        = A;
			this.b        = B;
			this.ra       = new Vector();
			this.rb       = new Vector();
			this.axis     = new Vector();
			this.normal   = new Vector();
			this.impulse  = new Vector();
			this.dv       = new Vector();
			this.va       = new Vector(va[0] * PHY2D.zoom, va[1] * PHY2D.zoom);
			this.vb       = new Vector(vb[0] * PHY2D.zoom, vb[1] * PHY2D.zoom);
			this.pa       = new Vector();
			this.pb       = new Vector();
			this.invDenom = 0;
			this.impulseN = 0;
			return this;
		}
		
		Constraint.prototype.presolve = function ()
		{
			// transform constraint points
			this.pa.transform(this.va, this.a.matrix);
			this.pb.transform(this.vb, this.b.matrix);
			
			// projection axis
			this.axis.sub(this.pb, this.pa);
			this.normal.scale(this.axis, 1 / (this.axis.len() + 0.0001));
			
			// calculate radius arms
			this.ra.sub(this.pa, this.a.matrix.pos).perp();
			this.rb.sub(this.pb, this.b.matrix.pos).perp();
			
			// compute denominator in impulse equation
			var ran = this.ra.dot(this.normal), rbn = this.rb.dot(this.normal);
			this.invDenom = 1 / (this.a.invMass + this.b.invMass + (ran * ran * this.a.invI) + (rbn * rbn * this.b.invI));
		}
		
		Constraint.prototype.solve = function ()
		{
			// calculate relative velocity in the axis, we want to remove this
			this.dv.sub(
				v1.addScale(this.b.vel, this.rb, this.b.angularVel),
				v2.addScale(this.a.vel, this.ra, this.a.angularVel)
			);
			
			// accumulated impulse
			var dist = this.axis.dot(this.normal);
			var newImpulse = (this.dv.dot(this.normal) + dist / kTimeStep) * this.invDenom + this.impulseN;
			
			// apply impulse
			this.impulse.scale(this.normal, newImpulse - this.impulseN)
			this.applyImpulse();
			this.impulseN = newImpulse;
		}
		
		Constraint.prototype.applyImpulse = Contact.prototype.applyImpulse;
		
		// manipulate constructor
		
		function Manipulate (B)
		{
			v1.sub(B.matrix.pos, pointer.pos);
			this.b        = B;
			this.rb       = new Vector();
			this.axis     = new Vector();
			this.normal   = new Vector();
			this.impulse  = new Vector();
			this.dv       = new Vector();
			this.vb       = new Vector().rotateIntoSpaceOf(v1, B.matrix);
			this.pb       = new Vector();
			return this;
		}
		
		Manipulate.prototype.solve = function ()
		{
			this.pb.transform(this.vb, this.b.matrix);
			this.axis.sub(this.pb, pointer.pos);
			this.normal.scale(this.axis, 1 / (this.axis.len() + 0.0001));
			this.rb.sub(this.pb, this.b.matrix.pos).perp();
			this.dv.addScale(this.b.vel, this.rb, this.b.angularVel);
			this.impulse.scale(this.normal, this.dv.dot(this.normal) + this.axis.dot(this.normal));
			this.b.vel.sub(this.b.vel, this.impulse);
			this.b.angularVel -= this.impulse.dot(this.rb) * this.b.invI;
		},

		// =============== Interface ==================
		
		this.render = function ()
		{
			// brute force AABB broadphase
			contacts.index = 0;
			for (var i = 0, len = objects.length; i < len - 1; i++)
			{
				var A = objects[i];
				for (var j = i + 1; j < len; j++)
				{
					var B = objects[j];
					if (A.invMass || B.invMass)
					{
						var a = A.aabb,
							b = B.aabb;
						if (
							Math.abs(b.x - a.x) - (a.w + b.w) < 0 &&
							Math.abs(b.y - a.y) - (a.h + b.h) < 0
						) A.contact(B);
					}
				}
			}


			
			// solver loop
			var len = contacts.index;
			var cln = joints.length;
			
			// solve contacts
				
			for (var i = 0; i < cln; i++)
			{
				joints[i].presolve();
			}
			
			for (var j = 0; j < 5; j++) // numIterations
			{
				// solve constraints
				for (var i = 0; i < cln; i++)
				{
					joints[i].solve();
				}
				// solve contacts
				for (var i = 0; i < len; i++)
				{
					contacts[i].solve();
				}
			}

			// integration loop
			for (var i = 0, len = objects.length; i < len; i++)
			{
				objects[i].integrate();
			}
				
			// draw
			for (var i = 0; i < len; i++)
			{
				var rb = objects[i];
				rb.draw();
			}
			
			// manipulate
			
			if (manipulate)
			{
				manipulate.solve();
				ctx.beginPath();
				ctx.lineWidth = 2;
				ctx.setLineDash([2,2]);
				ctx.moveTo(pointer.pos.x, pointer.pos.y);
				ctx.lineTo(manipulate.pb.x, manipulate.pb.y);
				ctx.strokeStyle = "rgb(255,255,255)";
				ctx.stroke();
				ctx.closePath();
				ctx.fillStyle = "rgb(255,255,255)";
				ctx.arc(manipulate.pb.x, manipulate.pb.y, 5, 0, 2 * Math.PI);
				ctx.fill();
			}
		}
		
		// create new rectangles

		this.rectangle = function (x, y, w, h, mass, angle, img)
		{
			x *= this.zoom;
			y *= this.zoom;
			w *= this.zoom;
			h *= this.zoom;
			var vertices = [
				w / 2, -h / 2, -w / 2, -h / 2, -w / 2, h / 2, w / 2, h / 2
			];
			return new Polygon(x, y, w, h, vertices, mass, angle, img);
		}
		
		this.triangle = function (x, y, w, h, mass, angle, img)
		{
			x *= this.zoom;
			y *= this.zoom;
			w *= this.zoom;
			h *= this.zoom;
			var vertices = [
				w / 2, h / 2, 0, - h / 2, -w / 2, h / 2
			];
			return new Polygon(x, y, w, h, vertices, mass, angle, img);
		}
		
		// add constraint
		
		this.addJoint = function (A, va, B, vb)
		{
			joints.push(
				new Constraint(A, va, B, vb)
			);
		}
		
		// delete static objects

		this.deleteStatic = function ()
		{
			var k = objects.length;
			while (k--)
			{
				var p = objects[k];
				if (!p.invMass)
				{
					objects.splice(k, 1);
				}
			}
		}
		
		// global zoom
		this.zoom = 1;
		this.stopManipulate = function () { manipulate = null; }

	}).call(
		PHY2D,
		ctx,
		pointer,
		7.8, // gravity
		1/30, // time step
		0.0, // friction
		true, // motionAABB
		"https://raw.githubusercontent.com/Bowsse/Bowsse.github.io/master/" // path images
	);

	////////////////

	screen.resize();

	
	
	var x4 = 150/(days/100);
	var xw = 200/(days/100);
	var xh = 200/(days/100);
	
    var irb = 0;
	

	var F = PHY2D.rectangle(300, 155, 200, 200, 1, 55, {src: "molang1.png", w:320, h:400,y:-14});
	var F = PHY2D.rectangle(300, 155, 200, 200, 1, 55, {src: "molang2.png", w:320, h:400,y:-14});
	var F = PHY2D.rectangle(300, 155, 200, 200, 1, 55, {src: "molang3.png", w:320, h:400,y:-14});
	var F = PHY2D.rectangle(300, 155, 200, 200, 1, 55, {src: "molang4.png", w:320, h:400,y:-14});



	var F = PHY2D.rectangle(300, 155, 200, 200, 1, 55, {src: "molang1.png", w:320, h:400,y:-14});
	var F = PHY2D.rectangle(300, 155, 200, 200, 1, 55, {src: "molang2.png", w:320, h:400,y:-14});
	var F = PHY2D.rectangle(300, 155, 200, 200, 1, 55, {src: "molang3.png", w:320, h:400,y:-14});
	var F = PHY2D.rectangle(300, 155, 200, 200, 1, 55, {src: "molang4.png", w:320, h:400,y:-14});


	var F = PHY2D.rectangle(300, 155, 200, 200, 1, 55, {src: "molang1.png", w:320, h:400,y:-14});
	var F = PHY2D.rectangle(300, 155, 200, 200, 1, 55, {src: "molang2.png", w:320, h:400,y:-14});
	var F = PHY2D.rectangle(300, 155, 200, 200, 1, 55, {src: "molang3.png", w:320, h:400,y:-14});
	var F = PHY2D.rectangle(300, 155, 200, 200, 1, 55, {src: "molang4.png", w:320, h:400,y:-14});



	var F = PHY2D.rectangle(300, 155, 200, 200, 1, 55, {src: "molang1.png", w:320, h:400,y:-14});
	var F = PHY2D.rectangle(300, 155, 200, 200, 1, 55, {src: "molang2.png", w:320, h:400,y:-14});
	var F = PHY2D.rectangle(300, 155, 200, 200, 1, 55, {src: "molang3.png", w:320, h:400,y:-14});
	var F = PHY2D.rectangle(300, 155, 200, 200, 1, 55, {src: "molang4.png", w:320, h:400,y:-14});





	while (irb < days)
	{
		var F = PHY2D.rectangle(500, 55, 51, 100, 1, 55, {src: "piu.png", w:100, h:100,y:-14});

		irb++;
	}
	
	
console.log(days);	
	// ==== main loop ====
	function run()
	{
		requestAnimationFrame(run);
		ctx.clearRect(0, 0, screen.width, screen.height);
		PHY2D.render();
	}
	// ==== start animation ====
	requestAnimationFrame(run);
}();
