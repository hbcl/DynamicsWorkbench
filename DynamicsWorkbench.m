(* ::Package:: *)

Off[General::spell1]


BeginPackage["DynamicsWorkbench`"];


(* ::Input:: *)
(*(* Version history:*)
(*v. 4.1 1/14: Altered ExportEOM to specify mFileName  *)
(*v. 4.0 7/13: improved Zee functions, added ExportEquations, ExportEOM, bug fix for gimbal jnt*)
(*v. 3.9 2/15/13: added drawing commands*)
(*v. 3.8 3/27/05: added ExportCode, ExportExpression *)
(*v. 3.7 3/15/04: added CenterDot notation, sped up slow dot products*)
(*v. 3.6 1/25/04: altered Norm to make compatible with Mma 5.0*)
(*v. 3.5 7/17/01: added facility to memorize dot products*)
(*v. 3.4 5/14/00: made slight change to Slider joint*)
(*v. 3.3 1/10/00: changed NewModel, added options*)
(* v. 3.2 10/2/99: Zee functions operational*)
(*v. 3.1 6/12/98: Zee functions in beta*)
(*v. 3.0 1/14/98: added Mma 3.0 functionality *)*)


(* ::Section:: *)
(*Help*)


AccCOM::usage = "AccCOM[body] returns a vector describing the acceleration of the center of mass of body with respect to the ground reference frame. AccCOM[{body1, body2,...}] returns the acceleration of the center of mass of the system consisting of the bodies listed.";


AccJnt::usage = "AccJnt[body] returns the velocity of the point on the inboard body which corresponds to the location of the joint connecting the two.";


AccJnt2::usage = "AccJnt2[body] returns the acceleration of the point on body which instaneously corresponds to the location of the joint between the body and the inboard body. Used in slider joints.";


AddFrame::usage = "AddFrame[frame,baseframe,jointtype,opts] adds a new reference frame to the frame tree. This frame is specified relative to a base frame which precedes it in the frame tree. Several Joint types are available.";


AddBody::usage = "AddBody[body,inboard,jointtype,opts] adds a new body to the system. This body is connected to the inboard body through a joint specified by jointtype. Several Joint types are available (type '?Joints'). Options InbToJnt and BodyToJnt specify vectors from the center of mass of the inboard body to the joint, and from the new body's center of mass to the joint. These vectors may be in any reference frame(s). Options Mass and Inertia describe physical parameters.  Other options specify the axes corresponding to the joint.";


AngAcc::usage = "AngAcc[body] returns the angular acceleration of the body with respect to ground.";


AngMom::usage = "AngMom[body,pnt] returns the angular momentum of the body about a given point, with respect to the ground reference frame. AngMom[body,pnt,ref] returns the angular momentum in a specified (inertial) frame. A list of bodies can also be provided. The given point pnt is treated as fixed in the reference frame. For angular momentum about a moving point, it is necessary to add the appropriate inertial term manually. For a single body, AngMom[body] returns the angular momentum about the body's center of mass.";


AngVel::usage = "AngVel[body] returns the angular velocity of the body with respect to ground. AngVel[body,frme] returns the angular velocity of the body with respect to frme";


AppFrc::usage = "AppFrc[body,force,point] applies vector force to
the body at a vector point from the center of mass.  Vectors may
be specified in defined reference frame.";


AppTrq::usage = "AppTrq[body,torque] applies a moment specified
by the vector torque to body.  The vector may be specified in
any defined reference frame.";


AtT0::usage = "AtT0[expr] returns expr set at time zero.";


Axis::usage = "Axis->axis specifies a rotational axis, and may be in the form of an axis number (1 through 3), a vector (e.g., ground[1]), or a list of unit vector coefficients (direction cosines). For joints with multiple rotatory axes, they are numbered Axis1, Axis2, etc.";


Ball::usage = "A Ball joint is a three degree-of-freedom rotatory joint. Its home position is coincident with the base reference frame (the frame of the inboard body). The generalized coordinates Qdof1 through Qdof4 describing the rotations are in the form of euler parameters {e1,e2,e3,e4}, where e4 is the cosine of the half-rotational angle.  The generalized speeds Udof1 through Udof3 are the angular velocities about the body's reference frame axes.";


Bodies::usage = "Bodies returns a list of the currently defined bodies.";


BodyToJnt::usage = "BodyToJnt[body] is the vector from the body's center of mass to the joint connecting it to the inboard joint.";


CastMtx::usage = "CastMtx[from,to] returns the transformation
matrix needed to cast a vector from one reference frame to
another.";


CastV::usage = "CastV[vector,frame] will cast any vector into
a specified reference frame.";


CollectE::usage = "CollectE[eom] puts the equations of motion into a form simplified with respect to the udots, so that they do not appear multiple times. CollectE[eom, x] simplifies with respect to any expression x appearing in the equations.";


ConvertList::usage = "ConvertList[vector, frame] is a utility function which converts either a vector, a list of unit vector coefficients (direction cosines), or the number between 1 and 3 designating an axis into the form of a list of the unit vector coefficients for the specified frame.";


Cross::usage = "Cross[vector1,vector2] or
(vector1 ~X~ vector2) performs the cross product of
two vectors.";


D::usage = "Partial derivatives of single or multiple vectors are found using D[ vector, scalar, in->frame, NonConstants->{functions} ], where in and NonConstants are options. Option in is used to specify the reference frame the partial derivative is desired in, and defaults to ground.";

Dt::usage = "Total derivatives of a single or multiple vector are found using Dt[ vector, scalar, in->frame, Constants->{constants} ], where in and Constants are options. Option in is used to specify the reference frame the total derivative is desired in, and defaults to ground. If the angular velocities of all of the reference frames vector is defined in are themselves defined, the total time derivative is computed using those angular velocities.";


Dist::usage = "Dist[point1,point2] returns the distance
between two points, which need not be described in the
same reference frames.";


Dot::usage = "Dot[vector1,vector2] or vector1 . vector2
performs the dot product of two vectors (or a vector
and a dyadic).";

Drawing::usage = "Drawing[body] or Drawing[{body1,...}] yields a list of graphic objects that connect the joints and show the center of mass positions.";


Dyadic::usage = "A dyadic may be performed by using
vector1 ** vector2.";


EOM::usage = "EOM[] generates the equations of motion. EOM[Simplify->On] applies simplification at low levels to attempt to generate more compact equations, but at the expense of more execution time.";


FindEquil::usage = "FindEquil[ eqn, {x,x0} ] finds the equilibrium point in the given equation, for variable x. The initial guess of x0 is used to start a Newton search. For a list of equations, use FindEquil[ eqns, {x,x0}, {y,y0},...] to find equilibrium points in several variables x, y,...";


Fixed::usage = "A Fixed joint is a zero degree-of-freedom joint. Its axes can be offset from those of the base reference frame using the Offset option. Offset may be a matrix of direction cosines or a list of vectors specifying the offset axes in terms of the base frame.  Examples: Offset->{{0,1,0},{1,0,0},{0,0,1}}, Offset->{body[2],body[1],body[3]}. A Fixed joint has no generalized speeds or coordinates associated with it.";


Force::usage = "Force[body] returns a list of the forces acting on the body.";


Frames::usage = "Frames returns a list of the currently defined reference frames.";


GActFrc::usage = "GActFrc[body,n] returns the generalized active force of body with respect to the nth generalized speed. GActFrc[body] returns the generalized active forces for body with respect to each of the generalized speeds. GActFrc[All] returns the sum of the generalized active forces for all bodies, with respect to each of the generalized speeds.";


Gimbal::usage = "A Gimbal joint is a three degree-of-freedom rotatory joint. Its home position is coincident with the base reference frame (the frame of the inboard body). The generalized coordinates describing the rotations are the successive rotation angles about the three specified rotation axes, described by options Axis1, Axis2, and Axis3. Use Ball if euler parameters are preferred as generalized coordinates. Axes may be specified by coordinate number, or a list of direction cosines. Examples: Axis1->3, Axis2->1, Axis3->2; Axis1->{0,0,1}, Axis2->{0,.707,.707}, Axis3->{.707,0,.707}.";


GInerFrc::usage = "GInerFrc[body,n] returns the generalized
inertia force of body with respect to the nth generalized
speed. GInerFrc[body] returns the generalized inertia forces
for body with respect to each of the generalized speeds.
GInerFrc[All] returns the sum of the generalized inertia
forces for all bodies, with respect to each of the
generalized speeds.";


ground::usage = "ground is the default Newtonian reference frame, which all other frames are defined relative to.";


Hinge::usage = "A Hinge joint is a single degree-of-freedom rotatory joint. The axis about which the joint rotates is specified by the option Axis, which may be the coordinate number, a list of direction cosines, or a vector in the base reference frame (the frame of the inboard body).  Examples: Axis->3, Axis->{0,0,1}, Axis->body[1]. The generalized speed, Udof, is defined to be the time-derivative of generalized coordinate Qdof.";


InbToJnt::usage = "InbToJnt[body] is the vector from the center of mass of the inboard body to the joint connecting it to the argument body.";


Inboard::usage = "Inboard[body] returns the name of the body inboard to the argument.";


Inertia::usage = "Inertia[body] is the inertia dyadic of the body.";


JntToJnt::usage = "JntToJnt[body] is the vector from the joint of the inboard body to the joint of the outboard body. This vector is used internally with slider joints.";


Joints::usage = "The following joint types are supported: Hinge, Ball, Slider, Fixed, Gimbal, SixDOF, UJoint. For more information, use ? followed by the joint type.";


Kids::usage = "Kids[frame] returns a list of reference frames below the argument's entry in the frame tree.";


Kinematics::usage = "Kinematics returns the set of kinematical differential equations.";


LinMom::usage = "LinMom[body] is the linear momentum of the body. \
LinMom[{body1,body2,...}] is the linear momentum of the system \
consisting of the bodies listed."; 


Mass::usage = "Mass[body] is the mass of the body.";


MassMatrix::usage = "MassMatrix[eom] returns the mass matrix and other terms from the specified equations of motion."; 


NewModel::usage = "NewModel[] initializes the Dynamics
Workbench for a new model. The generalized coordinates
and generalized speeds, and information regarding bodies
and reference frames are all cleared. NewModel[Zees->On] indicates that the subsequent model is to be constructed with intermediate terms stored in the 'zees' variable";


Nonholonomic::usage = "Nonholonomic returns a list of the nonholonomic constraints.  These include holonomic constraints rising from closed-loop mechanisms, which are treated as if they were nonholonomic and therefore expressed in differential form.";


Offset::usage = "Offset->axeslist is an option to AddFrame which offsets the new reference frame with respect to the base frame.  The default is no offset.  Offset is specified as the coordinates of each of the new frame's axes, in terms of the base frame.  These may specified as vectors, a list of three lists of direction cosines, or a list of three axis numbers.";


Parents::usage = "Parents[frame] returns a list of reference frames above the argument's entry in the frame tree.";


PosCOM::usage = "PosCOM[body] returns a vector describing the position of the center of mass of body with respect to the ground reference frame. PosCOM[{body1,body2,...}] returns the position of the center of mass of the system consisting of the bodies listed.";


PosPnt::usage = "PosPnt[point,body] returns the vector
describing the position of point (entered as a vector
relative to the center of mass) fixed to body.";


PrtVel::usage = "PrtVel[expr,n] or PrtVel[expr,u[n]] returns
the partial velocity of expression with respective to the
nth generalized speed.";


PV::usage = "PV[{c,frame,n}..}] is a head used by the Dynamics Workbench to denote vectors. Each 'Packaged Vector' is placed within this head and consists of a series of lists denoting a constant multiplied by a unit vector n of the reference frame: c frame[n]. PV's are not normally visible to users.";


q::usage = "q[n] refers to the nth generalized coordinate.";


qdofs::usage = "qdofs is a list of the currently defined generalized coordinates.";


qdots::usage="qdots[qs] returns a list of qdots, of the form q[1]', q[2]' etc. for the specified coordinates.";


RelativeTo::usage = "RelativeTo->frame is an option for AddBody which specifies an alternate reference frame which the new body's frame is referenced to. The default for RelativeTo is the inboard body's frame, so that the generalized coordinate q[] measures a displacement relative to the inboard body. In 2D systems, it is sometimes helpful to define RelativeTo->ground to yield an absolute angle in space. In 3D, there is not a similarly meaningful definition for absolute angles, and RelativeTo should not be used.";


ResltF::usage = "ResltF[body] returns the resultant of the forces
acting on the body.";


ResltT::usage = "ResltT[body] returns the resultant of the
torques or moments acting on the body.";


Rstar::usage = "Rstar[body] returns the inertia forces acting on
body.";


SixDOF::usage = "A SixDOF joint is a six degree-of- freedom joint which allows for full translational and rotational movement between bodies.  Its location ahd home position are specified by the BodyToJoint option. The translational axes are specified by TAxis1, TAxis2, and TAxis3. The translations along those axes are specified by generalized coordinates Qdof1, Qdof2, and Qdof3, and by generalized speeds Udof1, Udof2, and Udof3. The rotational orientation is described by four euler parameters, specified by generalized coordinates Qdof4 through Qdof7. The angular velocities about the base frame's axes are specified by generalized speeds Udof1, Udof2, and Udof3.";  


Slider::usage = "A Slider joint is a single degree-of- freedom prismatic joint. Its location and home position are specified by the BodyToJoint option. The translational axis is specified by TAxis, which may be a coordinate number, a list of direction cosines, or a vector in the base reference frame.  Examples: TAxis->3, TAxis->{0,0,1}, TAxis->body[1]. The generalized coordinate, Qdof, is defined as a position sliding relative to the inboard body. The generalized speed, Udof, is defined to be the time-derivative of generalized coordinate Qdof.";


SolveEqn::usage = "SolveEqn[eqn, vars] solves a list of equations in eqn for the given variables, and puts them in order as an == expression. SolveEqn is useful if the -> expression returned by Solve is inconvenient.";


SolveKinematics::usage="SolveKinematics solves the kinematical differential equations and returns the replacements of qdots with u's. This is useful for cases such as taking the derivative of an expression, yielding q[1]' and automatically converting it into u[1]."


States::usage = "States returns a list of all of the states
currently used in the model."; 


t::usage = "t is time."; 


TAxis::usage = "TAxis->taxis specifies a translational axis, and may be in the form of an axis number (1 through 3), a vector (e.g., ground[1]), or a list of unit vector coefficients (direction cosines). For joints with multiple translational axes, they are numbered TAxis1, TAxis2, etc.";


Tmtx::usage = "Tmtx[frame] is the transformation matrix which maps between frame and its inboard reference frame. Tmtx returns a list of coordinates of the three axes of frame, each in the form of a sublist comprising coefficients for axes of reference frame.";


Torque::usage = "Torque[body] returns a list of the torques acting on the body.";


Tstar::usage = "Tstar[body] returns the inertia torque acting on
body.";


u::usage = "u[n] refers to the nth generalized speed.";


udofs::usage = "udofs is a list of the currently defined generalized speeds.";


udots::usage="udots[us] returns a list of udots, of the form u[1]', u[2]', etc. for the specified speeds.";


UJoint::usage = "A Ujoint is a two degree-of-freedom rotatory joint. Its home position is coincident with the base reference frame (the frame of the inboard body). The generalized coordinates describing the rotations are the successive rotation angles about the two specified rotation axes, Axis1 and Axis2. Axis1 should be treated as being fixed in the inboard body frame, and Axis2 in the successive frame. It is recommended that the axes be specified in measure number form, e.g. Axis1->{0,0,1}, Axis2->{0,1,0} means successive z, y rotations starting in the inboard frame. Or equivalently, using axis numbers, e.g. Axis1->3,Axis2->2, meaning inboard frame axes 3 and then 2. The generalized speeds, Udof1 and Udof2, are defined to be the time-derivatives of the respective generalized coordinates, Qdof1 and Qdof2.";


VelCOM::usage = "VelCOM[body] returns a vector describing the velocity of the center of mass of body with respect to the ground reference frame.  VelCOM[{body1, body2,...}] returns the velocity of the center of mass of the system consisting of the bodies listed.";

VelPntFixed::usage = "VelPntFixed[comToPoint, body] returns the velocity of a point fixed on the body, where COMToPoint is the vector from COM to the point.";

AccPntFixed::usage = "AccPntFixed[comToPoint, body] returns the acceleration of a point fixed on the body, where COMToPoint is the vector from COM to the point.";


VelJnt::usage = "VelJnt[body] returns the velocity of the point on the inboard body which corresponds to the location of the joint connecting the two.";


VelJnt2::usage = "VelJnt2[body] returns the velocity of the point on body which instantaneously corresponds to the location of the joint between body and the inboard body. Used in slider joints, where VelJnt (from inboard body) and VelJnt2 differ by the sliding velocity.";


XDot::usage = "XDot[eom, state] returns the state derivative evaluated using the provided equations eom, and the current state.  The command States describes the ordering of generalized coordinates and speeds in the state vector.";


Z::usage = "Z[n] refers to the nth zee variable, but does
not return the value. Use zees[[n]] to retrieve the value
of the intermediate expression. Use NewModel[Zees->On] to specify that zees should be formed.";


Zees::usage="Zees->On is an option to NewModel that specifies that intermediate 'zee' variable should be generated in order to make expressions shorter.  Z[n] refers to the nth zee variable, but does not return the value. Use zees[[n]] to retrieve the value of the intermediate expression. Use NewModel[Zees->On] to specify that zees should be formed.";
Zee::usage="shortexpr = Zee[expr] finds factors that are linear in u or u', and replaces them with shorthand Z's. Returns a shorter version of the expression, also appends new shorthand items in zees. Similarly, finds constant terms (independent of t) and stores in zeesC.";
UnZee::usage="longexpr = UnZee[expr] restores the longhand form of an expression, removing the shorthand Z's.";
MinimumZee::usage="{nzeesC,nzees} = MinimumZee[expr] returns lists of the constant zees and regular zees indices that appear in the given expression or expressions. It leaves out the zees that do not appear.";


formtrigreplacements::usage="formtrigreplacements[equations] returns a list of shorthand trigonometric replacements that could be used to condense the equations given. The output is of the form {Sin[q[1]]->\"s1\",Cos[q[2]]->\"c2\"}. Options include ArgReplacements->{{gamma,\"g\",alpha,\"alp\"}} for shortening other argument variables, and TrigReplacements->{{Sin,\"s\",{Cos,\"c\"}} for changing the default trigonometric function substitutions.";


formqureplacements::usage=
"formqureplacements returns a list of shorthand replacments for the q's, u's, and udot's, e.g. {q[1]\[Rule]\"q1\",u[1]->\"u1\",u[1]'->\"u1dot\"}.";


ExportExpression::usage="ExportExpression[expression, \"name\"] formats expression as Matlab code, producing a string that begins with 'name = ' followed by the expression. ExportExpression can export a matrix, a symmetric matrix, a list, or a scalar. Formatting options include EndLine\[Rule]\"\\n\", Continuator\[Rule]\"...\", LineLength\[Rule]78, Deliminator\[Rule]\" \", which specify the end-of-line character, line continuator, maximum line length of the output, and the deliminators where line breaks can occur (default to space character). The options to formtrigreplacements can also be used.";


ExportCode::usage="\!\(\*
StyleBox[\"ExportCode\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\" \",\nFontWeight->\"Bold\"]\)writes a Matlab function to evaluate the state-derivative. Features include substitution of q's and u's for the state x, and substitution of shorthand trig functions such as s1 for sin(q(1)). ExportCode[{massmatrix,rhs}] uses the mass matrix and right-hand side given (usually found using MassMatrix[]). Alternatively, ExportCode[equations] will call MassMatrix[] automatically. Several optional arguments may be given. Forces\[Rule]{T1,T2} specifies forces that enter linearly into the equations; ExportCode will export a matrix of coefficients for these forces. Another option is Expressions\[Rule]{{KE,\"KE\"},{PE,\"PE\"}}, which outputs the specified additional expressions, named with the strings given. Each expression must be paired with a name within a list. The output of ExportCode can be sent to a file with OutputFile\[Rule]\"filename\", which is automatically produced in the current directory. Formatting options include EndLine\[Rule]\"\\n\", Continuator\[Rule]\"...\", LineLength\[Rule]78, Delimiter\[Rule]\" \", which specify the end-of-line character, line continuator, maximum line length of the output, and the delimiters where line breaks can occur (default to space character). The options to formtrigreplacements can also be used.";
ExportEOM::usage="\!\(\*
StyleBox[\"ExportEOM\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\" \",\nFontWeight->\"Bold\"]\)writes a Matlab function to evaluate the state-derivative. Features include substitution of q's and u's for the state x, and substitution of shorthand trig functions such as s1 for sin(q(1)). ExportEOM[{massmatrix,rrhs}] uses the mass matrix and right-hand side given (usually found using MassMatrix[]). Alternatively, ExportEOM[equations] will call MassMatrix[] automatically. Several optional arguments may be given. Forces\[Rule]{T1,T2} specifies forces that enter linearly into the equations; ExportEOM will export a matrix of coefficients for these forces. Another option is Expressions\[Rule]{{KE,\"KE\"},{PE,\"PE\"}}, which outputs the specified additional expressions, named with the strings given. Each expression must be paired with a name within a list. The output of ExportEOM can be sent to a file with OutputFile\[Rule]\"filename\", which is automatically produced in the current directory. The Matlab function is named after the OutputFile if specified, with default \"fstatederivative\". It can also be explicitly set with option mFileName->\"name\". Formatting options include EndLine\[Rule]\"\\n\", Continuator\[Rule]\"...\", LineLength\[Rule]78, Delimiter\[Rule]\" \", which specify the end-of-line character, line continuator, maximum line length of the output, and the delimiters where line breaks can occur (default to space character). The options to formtrigreplacements can also be used.";


LinearTerms::usage="LinearTerms[expressions,terms] returns a matrix of linear coefficients of the variables listed in terms,that appear in the expressions (also a list) given. One option is Simplifier\[Rule]Simplify, which causes the specified function to be applied to the matrix. Also see UnLinearTerms.";
UnLinearTerms::usage="UnLinearTerms[expressions,terms] returns all of the remaining terms that are not linear in the LinearTerms. One option is Simplifier->Simplify, which causes the specified function to be applied to the output.";


fmtString::usage="fmtString[string] reformats a long String argument and tries to break it into a bunch of lines with a Matlab (or other) continuator,so that the output is fairly readable code. Options include EndLine\[Rule]\"\\n\" (end of line character), Continuator\[Rule]\"...\" (Matlab or other line continuator), LineLength\[Rule]78 (max line length), Delimiter\[Rule]\" \" (look for spaces to break lines on; lists of delimiters are also acceptable here).";


SpatialVelC::usage="SpatialVelC[body] finds the 6D spatial velocity-c (angular velocity and COM velocity, in body-fixed coords) for a body. This is a 6-tuple in terms of q's and u's: {omega1, omega2, omega3, vc1, vc2, vc3}. Also can be called with a list of bodies, or SpatialVelC[Bodies] or SpatialVelC[All] for all the bodies in a system. See also SpatialVel for spatial velocity for inboard joint.";SpatialJacobianC::usage="SpatialJacobianC[body] finds the Jacobian of the body pose c, where the rows are the angular velocities and COM velocity (6-tuple) for a body (in body-fixed coords), and there is one column for each generalized speed u. It is equivalent to the partial derivative of the velocities with respect to the u's. Also can be called with a list of bodies, or SpatialJacobianC[Bodies] or SpatialJacobianC[All] for all the bodies in a system, with 6N rows for N bodies. See also SpatialJacobian for spatial jacobian for inboard joint.";
SpatialVel::usage="SpatialVel[body] finds the 6D spatial velocity (angular velocity and inboard joint velocity, in body-fixed coords) for a body. This is a 6-tuple in terms of q's and u's: {omega1, omega2, omega3, vo1, vo2, vo3}. Also can be called with a list of bodies, or SpatialVel[Bodies] or SpatialVel[All] for all the bodies in a system. See also SpatialVelC for spatial velocity for center of mass.";SpatialJacobian::usage="J=SpatialJacobian[body] finds the Jacobian of the body pose, where the rows are the angular velocities and inboard joint velocity (6-tuple) for a body (in body-fixed coords), and there is one column for each generalized speed u. It is equivalent to the partial derivative of the velocities with respect to the u's. Also can be called with a list of bodies, or SpatialJacobian[Bodies] or SpatialJacobian[All] for all the bodies in a system, with 6N rows for N bodies. See also SpatialJacobianC for spatial jacobian for center of mass.";
ConstraintMatrixC::usage="ConstraintMatrixC[body] returns a matrix that denotes the constraints enforced by a body's inboard body. The matrix has as many rows as the constraint entails, and as many columns as the number of spatial velocities (in terms of center of mass velocity, 6 times the total number of bodies). The constraints are enforced through (matrix)*v = 0, where v is a stack of spatial velocities. Also can be called with a list of bodies or All for all bodies. Optionally can also be called with ConstraintMatrixC[body, jointtype] to specify a particular type of joint. By default it uses the inboard joint used to add the named body.";
SpatialInertiaC::usage="SpatialInertiaC[body] returns the spatial inertia matrix for the body about its center of mass. The matrix is a block diagonal with blocks consisting of a 3x3 central inertia matrix and 3x3 translational mass matrix for each body. Also can be called with a list of bodies or All for all bodies.";
InertiaToMatrix::usage="InertiaToMatrix[body] converts an inertia dyadic into an inertia matrix with a given basis frame. InertiaToMatrix[body] returns a body's inertia dyadic in the body-fixed frame. InertiaToMatrix[inertiaDyadic,frame] converts an arbitrary dyadic into a matrix based on the given frame.";
SpatialVelCToUMatrix::usage="SpatialVelCToUMatrix[body] produces a matrix that converts spatial velocities for center of mass into the u's. The number of rows equals the number of bodies, and the number of columns is 6 times the total number of bodies in the system. Also works with a list of bodies, or All for all bodies. Also call SpatialVelCToUMatrix[body, jointtype] to construct the matrix for a particular joint type. This is a helper function used internally.";    
BodyTransform::usage="BodyTransform[vector,body] transforms a Dynamics Workbench vector into a body's coordinate frame and expresses it as the measure numbers in that frame. In other words, if it returns the list {x,y,z} for a body, those are the xyz components of the vector in the body's frame. The vector is equivalent to (x body[1] + y body[2] + z body[3]). Also use GroundTransform to automatically transform into ground/Newtonian frame.";
GroundTransform::usage="GroundTransform[vector] transforms a Dynamics Workbench vector into the ground frame and expresses it as the measure numbers of ground coords. In other words, if it returns the list {x,y,z} for a body, those are the xyz components of the vector in the ground frame. Also see BodyTransform.";


JntAxes::usage="JntAxes[body] returns the joint axes associated with the body. Returns a Dynamics Workbench vector. In most cases, a joint axis is defined in terms of the inboard body's frame. For measure numbers in the body frame, use JntAxesBodyFrame.";
JntAxesBodyFrame::usage="JntAxesBodyFrame[body] returns the joint axes associated with the body, as measure numbers in that body's frame. Equivalent to JntAxes[body] except returns a list of measure numbers instead of a Dynamics Workbench vector.";
JntType::usage= "JntType[body] returns the type of joint (e.g. Hinge, Fixed, Slider) associated with the body. Also can be applied to a list, e.g. JntType[Bodies].";


Crossify::usage="Skew-symmetric matrix repesentation of a cross product operator. Crossify[{x,y,z}] takes the coordinates of a vector v = {x,y,z} and returns the linear operator equivalent to 'v cross...' Also works with Dynamics Workbench vectors, which are automatically transformed into ground frame. For example, Crossify[ground[3]]. Or transformed into specified body frame, Crossify[ground[3],body].";


Clear[Axis1,Axis2,Axis3,TAxis1,TAxis2,TAxis3,
	Udof,Udof1,Udof2,Udof3,Udof4,Udof5,Udof6,
	Qdof,Qdof1,Qdof2,Qdof3,Qdof4,Qdof5,Qdof6,Qdof7
	];


(* ::Input:: *)
(*Begin["`Private`"];*)


(* ::Section:: *)
(*Simple Utility Functions*)


(* ::Special1:: *)
(*AngVel[ frme1, frme2 ] returns the angular velocity of frme1 relative to frme2*)


AngVel[ frme1_, frme2_ ] := AngVel[ frme1 ] - AngVel[ frme2 ];


(* ::Special1:: *)
(*Common returns the leftmost identical elements between two lists. This is used in traversing trees.*)


Common[a_List, b_List] := Module[
	{n = Min[ Length[a], Length[b] ]},
	Part[ a, Flatten[ Position[
		MapThread[ SameQ, {Take[a,n], Take[b,n]} ], True
		] ] ]
]


(* ::Special1:: *)
(*Uncommon1 takes list a, and removes the leftmost common elements between a and b.  Example: a = {x, y, z}, b = {x, t, u}, Uncommon1[ a, b ] = {y, z}. This is used in traversing trees.*)


Uncommon1[a_List, b_List] := Module[
	{n = Min[ Length[a], Length[b] ]},
	Delete[ a, Position[
		MapThread[ SameQ, {Take[a,n], Take[b,n]} ], True
		] ]
]


(* ::Special1:: *)
(*Uncommon2 takes list b, and removes the leftmost common elements between a and b.  Example: a = {x, y, z}, b = {x, t, u}, Uncommon1[ a, b ] = {t, u}*)


Uncommon2[a_List, b_List] := Module[
	{n = Min[ Length[a], Length[b] ]},
	Delete[ b, Position[
		MapThread[ SameQ, {Take[a,n], Take[b,n]} ], True
	] ]
]


(* ::Special1:: *)
(*Wrap is a hack that forces Map to work on an empty list, without which Map does not apply any function*)


Wrap[ {} ] = {{}};
Wrap[ a_List ] := a;


(* ::Special1:: *)
(*Mod3 performs modulus on the ring {1,2,3}*)


Mod3[ n_Integer ] := Mod[ n-1, 3] + 1;


(* ::Special1:: *)
(*Normed gives the normalized (magnitude 1) vector v*)


Normed[ v_ ] := v / PowerExpand[ Sqrt[ v.v ] ];


NormV[ v_ ] := PowerExpand[ Sqrt[ v.v ] ];


Normed[ x:PV[__] ] := x / NormV[ x ];


(* ::Text:: *)
(*Following is deprecated, old tester version of NormV.*)


(* ::Input:: *)
(*NormV[ x:PV[__] ] := PowerExpand[ Sqrt[ Simplify[ Apply[ Plus, Map[ #^2 &, *)
(*	Map[ #[[1]]&, CastV[x, ground ] ] ] ] ] ] ]; (* deprecated, make sure cell is not initialization, not evaluatable *)*)


(* ::Special1:: *)
(*CoordQ determines whether the input is between 1 and 3*)


CoordQ[ n_ ] := IntegerQ[n] && Abs[n] >= 1 && Abs[n] <=3;


(* ::Special1:: *)
(*TripleQ determines whether the input is a vector of length 3*)


TripleQ[ v_ ] := VectorQ[v] && Length[v] == 3;


(* ::Special1:: *)
(*PrincipalAxisQ determines whether the input is a vector of length 3 which is either {0,0,1}, {0,1,0}, or {1,0,0}*)


PrincipalAxisQ[ x:{_,_,_} ] := Count[x,1] == 1 &&
	Count[x,0] == 2;


(* ::Special1:: *)
(*VecToList and ListToVecconvert between list and packaged vector *)
(*representations of vectors. VecToList[ v, frame ] produces a list of the measure numbers of vector v, in the coordinates of frame. ListToVec[ {x,y,z}, frame ] produces a Dynamics Workbench vector in the specified frame.*)


VecToList[ vec_, frame_ ] := {vec . frame[1], vec . frame[2],
  vec . frame[3]};


ListToVec[ lis_, frame_ ] := lis[[1]] frame[1] +
  lis[[2]] frame[2] + lis[[3]] frame[3];


(* ::Special1:: *)
(*ConvertVec is used to convert any of the possible three forms for*)
(*vectors into the standard packaged vector form. ConvertVec[ 3, frame ] just means frame[3]. ConvertVec[ {0,0,1}, frame ] converts the list of measure numbers into a Dynamics Workbench packaged vector. ConvertVec can also operate on a list of measure number lists, e.g. ConvertVec[ {{0,0,1},{0,1,0}}, frame ].*)


ConvertVec[ coord_?CoordQ, frme_ ] := ListToVec[ 
	Sign[coord]*RotateRight[{0,0,1},Abs[coord]], frme ];
ConvertVec[ axis_?TripleQ, frme_ ] := 
	ListToVec[ Normed[axis], frme ];
ConvertVec[ x:PV[__], frme_ ] := 
	x 1/Sqrt[(x.frme[1])^2+(x.frme[2])^2+(x.frme[3])^2];
ConvertVec[l:{(_?TripleQ)..}, frme_ ] :=
	Map[ ConvertVec[#,frme]&, l ]


(* ::Special1:: *)
(*ConvertList is used to convert any of the possible three forms for*)
(*vectors into a list form in the specified frame, returning a unit vector (as are needed to specify axes of rotation). ConvertList[ 3, frame ] produces the measure numbers {0, 0, 1}. ConvertList[ {0,0,1}, frame ] also produces the measure numbers, which are the same as the input list except normalized. ConvertList[ vector, frame ] takes a Dynamics Workbench vector and produces measure numbers for the specified frame, again normalizing to yield a unit vector.*)


ConvertList[ coord_?CoordQ, frme_ ] :=
	Sign[coord]*RotateRight[{0,0,1},Abs[coord]];
ConvertList[ axis_?TripleQ, frme_ ] :=
	Normed[axis];
ConvertList[ x:PV[__], frme_ ] :=
	VecToList[
		x 1/Sqrt[(x.frme[1])^2+(x.frme[2])^2+(x.frme[3])^2], frme];


(* ::Special1:: *)
(*ConvertOffset is used to convert one of the two possible forms for a reference frame offset. Input can be a transformation matrix, or a list describing the coordinates of each of the offset axes in the base reference frame.*)


ConvertOffset[M_?MatrixQ,frme_] := M;
ConvertOffset[{x_?CoordQ,y_?CoordQ,z_?CoordQ},frme_]:={ConvertList[x,frme],ConvertList[y,frme],ConvertList[z,frme]};
ConvertOffset[{x:PV[__],y:PV[__],z:PV[__]},frme_] :=
	{VecToList[x,frme],VecToList[y,frme],VecToList[z,frme]};


(* ::Special1:: *)
(*ConvertInertia is used to convert any of the possible three forms for inertia matrices into a dyadic form in the specified frame*)


ConvertInertia[I_?VectorQ,frme_] := I[[1]] frme[1]**frme[1] +
	I[[2]] frme[2]**frme[2] + I[[3]] frme[3]**frme[3];
ConvertInertia[I_?MatrixQ,frme_] := Apply[ Plus, Flatten[ Outer[ 
	NonCommutativeMultiply, Array[frme,3], Array[frme,3] ] ] *
	Flatten[ I ] ];
ConvertInertia[x:PV[__],frme_] := x;
ConvertInertia[0|0.,frme_]:=PV[{PV[{0,frme,1}],frme,1}];


(* ::Special1:: *)
(*If no dof number is assigned, PickDof is used to assign one automatically, based on the lowest number which is not already taken*)


PickDof[dofs_List] :=
	First[ Complement[ Range[ Length[dofs]+1 ], dofs ] ];
PickDof[dofs_List, num_?NumberQ] :=
	Take[ Complement[ Range[ Length[dofs]+num ], dofs ], num ];


(* ::Special1:: *)
(*ReplaceDofs takes a list of new dofs, and assigns the ones marked Automatic with new dofs which do not appear in the useddofs or newdofs lists.*)


ReplaceDofs[ newdofs_, useddofs_ ] := replDofs[ newdofs, 
	Flatten[ Position[ newdofs, Automatic ] ],
	PickDof[ Union[ Complement[ newdofs, {Automatic}],useddofs], 
	  Count[ newdofs, Automatic ] ] ];


replDofs[ newdofs_, {}, {}] := newdofs;


replDofs[ newdofs_List, loc:{_}, val:{_}] :=
	ReplacePart[ newdofs, val[[1]], loc];


replDofs[ newdofs_, loc_List, val_List ] :=
	replDofs[ ReplacePart[ newdofs, val[[1]], loc[[1]] ],
		Drop[loc, 1], Drop[val, 1] ];


replacesin = {
	c_. Sin[x_] Cos[y_] + c_. Cos[x_] Sin[y_] -> c Sin[x+y],
	c_. Sin[x_] Cos[y_] + d_. Cos[x_] Sin[y_] /; d===-c-> c Sin[x-y]};


replacecos = {
	c_. Cos[x_] Cos[y_] + c_. Sin[x_] Sin[y_] -> c Cos[x-y],
	c_. Cos[x_] Cos[y_] + d_. Sin[x_] Sin[y_] /; d===-c-> c Cos[x+y]};


replace2SC = {a_. Cos[x_] + b_. Cos[x_] -> (a+b) Cos[x],
a_. Sin[x_] + b_. Sin[x_] -> (a+b) Sin[x]};


replaceSC2 = {a_. Cos[x_]^2 + a_. Sin[x_]^2-> a};


replaceUs = {c_. u[n_] + d_. u[n_] -> (c+d) u[n]};


(* ::Special1:: *)
(*AtT0 is a simple function to enable use of this*)


AtT0[ expr_ ] := expr /. t->0;
AtT0[ expr__ ] := Flatten[{expr}] /. t->0;


(* ::Special1:: *)
(*FindEquil finds an equilibrium point for one or more equations*)


FindEquil[ eqn_, args:{_, _}.. ] :=
	FindRoot[ Evaluate[eqn], args ]


(* ::Text:: *)
(*Mass returns the masses of a bunch of bodies*)


Mass[(bodies_)?ListQ] := Mass /@ bodies; 


(* ::Text:: *)
(*VelCOM and AccCOM can return the velocity of the center of mass of a collection of bodies*)


VelCOM[(bodies_)?ListQ] := 
  Simplify[Plus @@ (Mass /@ bodies*VelCOM /@ bodies)/
    Plus @@ Mass /@ bodies]


AccCOM[(bodies_)?ListQ] := Simplify[Plus@@(Mass/@ bodies*AccCOM /@ bodies)/ Plus @@ Mass /@ bodies]


(* ::Subsubsection:: *)
(*Optional common stuff*)


(* ::Input:: *)
(*(* This is another way to do the common functions above,*)
(*   but it makes fewer assumptions of the form of the*)
(*   input lists *)*)


(* ::Input:: *)
(*common[a_List, b_List] := Select[ a, MemberQ[b, #]& ]*)


(* ::Input:: *)
(*uncommon1[a_List, b_List] := Select[ a, FreeQ[b, #]& ]*)


(* ::Input:: *)
(*uncommon2[a_List, b_List] := Select[ b, FreeQ[a, #]& ]*)


(* ::Text:: *)
(*SymmetricQ tests if a matrix is symmetric.*)


(* ::Input:: *)
(*SymmetricQ[(m_List)?MatrixQ] := m === Transpose[m]*)


(* ::Section:: *)
(*Vectors*)


(* ::Special1:: *)
(*All vectors are packaged together and tagged with PV.*)


(*PV/: Times[c_ /; FreeQ[c,PV],PV[{a_,frme_,n_}] ] :=
	PV[{a c,frme,n}]; *)
Timespkg[c_,{x_,y_,z_}] := {c x, y, z};
PV/: Times[c_ /; FreeQ[c,PV], PV[args__] ] :=
	Map[ Timespkg[c,#]&, PV[args] ];
PV/: Plus[ PV[x:{_,_,_}..], PV[y:{_,_,_}..] ] := 
	PV[x,y];
PV[{0,_,_}] := 0;
PV[{0,_,_},rest__] := PV[rest];
PV[rest1__,{0,_,_},rest2__] := PV[rest1,rest2];
PV[r1___,{c1_ /; FreeQ[c1,PV],frme_,n_},r2___,
	 {c2_ /; FreeQ[c2,PV],frme_,n_},r3___] :=
	 PV[r1,r2,r3,{c1+c2,frme,n}];
(*PV[r1:{_,_,_}...,{c1_ /; FreeQ[c1,PV],frme_,n_},r2:{_,_,_}...,
     {c2_ /; FreeQ[c2,PV],frme_,n_},rest___] :=
	PV[{c1+c2,frme,n},r1,r2,rest] ;*)
SetAttributes[PV, Orderless];


(* ::Special1:: *)
(*Though vectors are tagged with a head PV, meaning that they are packaged, the output form is set so that vectors look the same way as they are input*)


If[ TrueQ[ $VersionNumber >= 3.],
	Format[outV[{c_,v_,n_}]] := 
		HoldForm[c Subscript[UnderBar[v],n] ];
	Format[outV[{1,v_,n_}] ]  := Subscript[UnderBar[v],n],
(*else *)
	Format[outV[{c_,v_,n_}] ] :=
		HoldForm[c Subscripted[v[n]] ];
	Format[outV[{1,v_,n_}] ]  := HoldForm[Subscripted[v[n]] ] ];
(*Format[PV[x:{0,_,_}]] := 0;*)
Format[PV[x:{_,_,_}..] ] := Apply[ Plus, Map[ outV, {x}] ];



(* ::Special1:: *)
(*CollectV does garbage collection of packaged vectors, combining coefficients of identical elements*)


SetAttributes[CollectV,Orderless]


CollectV[ PV[r1___,{c1_,frme_,n_},r2___,
	{c2_,frme_,n_},r3___] ]:=
  CollectV[ PV[r1,r2,r3,{c1+c2,frme,n}] ];


CollectV[ PV[x__] ] := PV[x];


(* ::Section:: *)
(*Casting into Reference Frames*)


(* ::Special1:: *)
(*CastMtx returns the transformation matrix used to cast from one frame to another. It goes up the branches of the frame tree from the from frame up to the common parent to the two frames, and then down the branches to the to frame.*)


CastMtx[ from_, to_ ] :=
 	Apply[ Dot,
		Map[ Tmtx, 
			Wrap[Reverse[Uncommon1[ Parents[to], Parents[from] ]] ] 
		]
	] ~ Dot ~
		Transpose[ Apply[ Dot,
		Map[ Tmtx, 
			Wrap[Reverse[Uncommon2[ Parents[to], Parents[from] ]] ]
		]
	] ] /. Join[replacesin,replacecos,replaceSC2,replace2SC]


(* ::Special1:: *)
(*CastV will cast any vector into a target reference frame. It does this by casting each element one at a time, and adding the results. Cast1V takes one element at a time, casts it into three elements of the target frame, and assembles the results in a PV*)


CastV[ PV[x__], trgtfrme_ /; MemberQ[Frames, trgtfrme] ]:= 
	Simplify[ Apply[ Plus,
		Map[ Cast1V[#,trgtfrme]&, {x}, 1] ] /. replaceSC2 ];
		
CastV[ 0, trgtfrme_ ] := 0;


Cast1V[ {c_,frme_,n_}, trgtfrme_ ] := Apply[PV,
	MapThread[List,{c*CastMtx[trgtfrme,frme][[n]],
		{trgtfrme,trgtfrme,trgtfrme},{1,2,3}}] ] ;


(* ::Text:: *)
(*CML is the castmatrixlist, which is a way to memorize all the transformation matrices as they are computed, so nothing needs to be repeated.*)


CML[from_,to_]:=Module[{frame1,frame2},{frame1,frame2}=Sort[{from,to}];
If[OrderedQ[{from,to}],Identity,Transpose][CML[frame1,frame2]=CastMtx[frame1,frame2]]]


(* ::Section:: *)
(*Dot Products*)


(* ::Subsection:: *)
(*Dot products using CenterDot, defined at the top level*)


(* ::Text:: *)
(*Definition of zero dot products, to make sure that vectors with length 0 or 0. evaluate to zero dot product. *)


(* ::Input:: *)
(*PV/:CenterDot[0|0.,PV[__]]=0;*)
(*PV/:CenterDot[PV[__],0|0.]=0;*)


(* ::Text:: *)
(*Vectors dotted with a dyadic, and vice versa*)


(* ::Input:: *)
(*PV/:CenterDot[x:PV[{_,_,_}],PV[{y:PV[{_,_,_}],frame2_,coord2_}]]:=PV[{CenterDot[x,y],frame2,coord2}];*)


(* ::Input:: *)
(*PV/:CenterDot[PV[{PV[{c1_,frame1_,coord1_}],frame2_,coord2_}],x:PV[{_,_,_}]]:=PV[{CenterDot[PV[{c1,frame2,coord2}],x],frame1,coord1}];*)


(* ::Text:: *)
(*Distributive property of dot products*)


(* ::Input:: *)
(*PV/:CenterDot[PV[x:{_,_,_}..],PV[y:{_,_,_}..]]:=Apply[Plus,Flatten[Outer[CenterDot,Map[PV,{x}],Map[PV,{y}]]]]/;Length[{x}]+Length[{y}]>2*)


(* ::Text:: *)
(*Dot products of two vectors on the same unit vector is just their product*)


(* ::Input:: *)
(*PV/: CenterDot[PV[{c1_,frame_,coord_}],PV[{c2_,frame_,coord_}]]:=c1*c2*)


(* ::Text:: *)
(*Dot product of two vectors in the same reference frame, but different axes, is zero*)


(* ::Input:: *)
(*PV/: CenterDot[PV[{c1_,frame_,coord1_}],PV[{c2_,frame_,coord2_}]]:=0*)


(* ::Special1:: *)
(*Dot products are performed by overloading the Dot operator for PV objects. The overloading implements a distributed function on DotV, which extracts the dot product of two elements from the transformation matrix*)


(* ::Subsection:: *)
(*CenterDot defined to refer to DotV*)


PV/:CenterDot[PV[x:{_,_,_}..],PV[y:{_,_,_}..]]:=Apply[Plus,Distribute[DotV[{x},{y}],List]]/.Join[replacesin,replacecos,replaceSC2,replace2SC];


PV/: CenterDot[ 0 | 0., PV[__] ] = 0;
PV/: CenterDot[ PV[__], 0 | 0.] = 0;
CenterDot[ 0 | 0., 0 | 0.] = 0;


(* ::Subsection:: *)
(*Dot overloaded to refer to DotV*)


Unprotect[ Dot ];
Dot[ 0 | 0., 0 | 0.] = 0;
Protect[ Dot ];


PV/: Dot[ PV[x:{_,_,_}..], PV[y:{_,_,_}..] ] :=
	Apply[Plus, Distribute[ DotV[{x},{y}], List ] ] /.
	Join[replacesin,replacecos,replaceSC2,replace2SC];


PV/: Dot[ 0 | 0., PV[__] ] = 0;
PV/: Dot[ PV[__], 0 | 0.] = 0;


DotV[x:{_,_,_},{PV[y:{_,_,_}],frame2_,coord2_}] :=
	PV[ {DotV[x,y], frame2, coord2} ];


(*DotV[{PV[y:{_,_,_}],frame2_,coord2_},x:{_,_,_}] :=
	PV[ {DotV[x,y], frame2, coord2} ]; *)
DotV[{PV[y:{c1_,frame1_,coord1_}],frame2_,coord2_},x:{_,_,_}]:=
    PV[{DotV[{c1,frame2,coord2},x],frame1,coord1}];


DotV[{c1_,frame1_?(MemberQ[Frames,#]&),coord1_},
     {c2_,frame2_?(MemberQ[Frames,#]&),coord2_}] :=
	   c1 c2 CML[frame1,frame2][[coord2,coord1]]


(* ::Section:: *)
(*Cross Products*)


(* ::Special1:: *)
(*Cross products are given in the reference frame of the second argument. This is performed by breaking the 2nd arg into elementary PVs, and using CrossV, which requires that its inputs are in the same frame. CrossV takes each element of the 1st arg and crosses with the elemental 2nd arg*)


If[ TrueQ[ $VersionNumber >= 3. ],
	Unprotect[Cross] ];
Cross[ PV[x:{_,_,_}..], PV[y:{_,_,_}..] ] :=
	Module[{
		lftarg = Map[ CastV[ PV[x], #]&, Transpose[{y}][[2]] ],
		rgtarg = Map[ PV, {y} ]},
		Apply[Plus, MapThread[ CrossV, {lftarg, rgtarg} ] ]
	];
	
Cross[ expr_, 0 | 0.] := 0;
Cross[ 0 | 0., expr_ ] := 0;
Cross[ 0 | 0., 0 | 0.] := 0;


CrossV[PV[{c1_,frme_,n1_}],
	PV[{c2_,frme_,n2_}]] :=
	PV[{{ c1 c2, frme, Mod3[n1+2]},
	 {-c1 c2, frme, Mod3[n1+1]},
	 {     0, frme, n1 }}[[Mod3[n2-n1]]] ] /.
	Join[replacesin,replacecos,replaceSC2,replace2SC];


CrossV[PV[x:{_,frme_,_}..],PV[y:{_,frme_,_}] ] :=
	Apply[ Plus, Map[ CrossV[ PV[#], PV[y] ]&, {x}] ];


CrossV[ 0 | 0., PV[__] ] := 0;
CrossV[ PV[__], 0 | 0. ] := 0;


X[x_,y_] := Cross[x,y];


(* ::Section:: *)
(*Dyadics*)


(* ::Special1:: *)
(*To assemble dyadics, it is necessary to use a NonCommutativeMultiply to hold two vectors in the right order: r=a[1]**a[2] can be used so that later, r.b[1] gives the right result*)


PV/: PV[{x1_,x2_,x3_}] ** PV[{y1_,y2_,y3_}] :=
	PV[{y1 PV[{x1,x2,x3}],y2,y3}]

PV/: PV[x:{_,_,_}..]**PV[y:{_,_,_}..] := Apply[ Plus,
	Distribute[Map[ PV, {x} ] ** Map[ PV, {y} ], List ] ]
	
PV/: PV[__] ** (0 | 0.) := 0;
PV/: (0 | 0.) ** PV[__] := 0;


(* ::Section:: *)
(*Spatial vectors (after Featherstone)*)


BodyTransform[v_,body_]:={v.body[1],v.body[2],v.body[3]};
GroundTransform[v_]:=BodyTransform[v,ground];
SpatialVelC[body_ /; MemberQ[Bodies,body]]:=Join[BodyTransform[AngVel[#],#],BodyTransform[VelCOM[#],#]]& [body]//Simplify;
SpatialVelC[All]:=SpatialVelC[Bodies];
SpatialVelC[bodies_List]:=Apply[Join,Map[SpatialVelC[#]&,bodies]];
SpatialJacobianC[bodies_]:=LinearTerms[SpatialVelC[bodies],u[#]&/@udofs];
SpatialJacobianC[All]:=SpatialJacobianC[Bodies];


SpatialVel[body_ /; MemberQ[Bodies,body]]:=Join[BodyTransform[AngVel[#],#],BodyTransform[VelJnt[#],#]]& [body]//Simplify;
SpatialVel[All]:=SpatialVel[Bodies];
SpatialVel[bodies_List]:=Apply[Join,Map[SpatialVel[#]&,bodies]];
SpatialJacobian[bodies_]:=LinearTerms[SpatialVel[bodies],u[#]&/@udofs];
SpatialJacobian[All]:=SpatialJacobian[Bodies];


(* ::Text:: *)
(*Produce skew-symmetric matrix representation of a cross product. This is overloaded so that a DW vector also works, where it is transformed into ground (or specified body) coordinates and then crossified.*)


Crossify[{x_,y_,z_}]:={{0,-z,y},{z,0,-x},{-y,x,0}};
Crossify[v:PV[__]]:=Crossify[GroundTransform[v]];
Crossify[v:PV[__],body_]:=Crossify[BodyTransform[v,body]];


JntType[All]:=JntType/@Bodies;
JntType[bodies_List]:=JntType/@bodies;
JntAxes[bodies_List]:=JntAxes/@bodies;
JntAxesBodyFrame[body_]:=BodyTransform[JntAxes[body],body];


SpatialInertiaC[bodies_List]:=SparseArray[Band[{1,1}]->{##}]& @@(SpatialInertiaC[#]&/@bodies); (* found on MMa stackexchange *)
SpatialInertiaC[All]:=SpatialInertiaC[Bodies];
SpatialInertiaC[body_]:=ArrayFlatten[{{InertiaToMatrix[body],0},{0,Mass[body]IdentityMatrix[3]}}];


InertiaToMatrix[inertiaDyadic_,frme_]:=Table[frme[i].inertiaDyadic.frme[j],{i,3},{j,3}];
InertiaToMatrix[frme_]:=Table[frme[i].Inertia[frme].frme[j],{i,3},{j,3}];


SpatialVelCToUMatrix[All]:=SpatialVelCToUMatrix[Bodies];
SpatialVelCToUMatrix[bodies_List]:=Apply[Join,SpatialVelCToUMatrix[#,JntType[#]]& /@bodies];
SpatialVelCToUMatrix[body_,Hinge]:=Module[{convMat,inbTrans,jntaxis},
convMat=ConstantArray[0,{1, 6Length[Bodies]}]; (* Hinge is 1 dof *)
inbTrans = CastMtx[Inboard[body],body];
jntaxis=BodyTransform[JntAxes[body],body];
convMat[[1,SpatialAngVelIndex[Inboard[body]]]]=-jntaxis.inbTrans;
convMat[[1,SpatialAngVelIndex[body]]]=jntaxis;
convMat];
SpatialVelCToUMatrix[body_,Fixed]:={};
SpatialVelCToUMatrix[body_,Slider]:=Module[{convMat,inbTrans,jntaxis},
convMat=ConstantArray[0,{1, 6Length[Bodies]}]; (* Slider is 1 dof *)
inbTrans = CastMtx[Inboard[body],body];
jntaxis=BodyTransform[JntAxes[body],body];
convMat[[1,SpatialVelIndex[Inboard[body]]]]=+jntaxis.inbTrans;
convMat[[1,SpatialVelIndex[body]]]=jntaxis;
convMat];
SpatialVelCToUMatrix[body_,Gimbal|Ball|UJoint|SixDOF]:=Module[{},Print["This function has not been implemented yet!"];{}];
SpatialIndex[ground]={};
SpatialAngVelIndex[ground]={};
SpatialVelIndex[ground]={};
SpatialIndex[body_]:=Flatten[6(Position[Bodies,body]-1)][[1]];
SpatialAngVelIndex[body_]:=SpatialIndex[body]+{1,2,3};
SpatialVelIndex[body_]:=SpatialIndex[body]+{4,5,6};


ConstraintMatrixC[All]:=ConstraintMatrixC[Bodies];
ConstraintMatrixC[bodies_List]:=Apply[Join,ConstraintMatrixC[#,JntType[#]]& /@bodies];
ConstraintMatrixC[body_]:=ConstraintMatrixC[body,JntType[body]];
(* Fixed joint constraint: *)ConstraintMatrixC[body_,Fixed]:=Module[{inbAngVelColumns,inbVelColumns,bodyAngVelColumns,bodyVelColumns,inbTrans,consMat},
(* cols in constraint matrix corresponding to the inboard body: *)inbAngVelColumns=(* where ground is not a body *)Flatten[If[Inboard[body]===ground,{{},{},{}}, (* else *)6(Position[Bodies,Inboard[body]][[1]][[1]]-1)]+{1,2,3}];
inbVelColumns=inbAngVelColumns+3; (* translational velocities are right after ang vels *)
bodyAngVelColumns=6(Position[Bodies,body][[1]][[1]]-1)+{1,2,3};
bodyVelColumns=bodyAngVelColumns+3;
inbTrans = CastMtx[Inboard[body],body];(* transform from inboard to body *)
consMat=ConstantArray[0,{6,6Length[Bodies]}]; (* 6 x N array *)
(* Constrain fixed joint translation: *)
       consMat[[{1,2,3},inbVelColumns]]=inbTrans; (* velocity of inboard COM *)
       consMat[[{1,2,3},inbAngVelColumns]]=inbTrans.Crossify[BodyTransform[-InbToJnt[body],Inboard[body]]]; (* omega cross r, for inboard COM to joint *)
        consMat[[{1,2,3},bodyVelColumns]]=-IdentityMatrix[3];(* -velocity of body COM *)consMat[[{1,2,3},bodyAngVelColumns]]=-Crossify[BodyTransform[-BodyToJnt[body],body]]; (* - omega cross r, for body COM to joint *)
(* Next constrain slider joint rotation (none allowed) *)
consMat[[{4,5,6},inbAngVelColumns]]=-inbTrans; (* (omega1 - omega2) = 0 *)
consMat[[{4,5,6},bodyAngVelColumns]]=IdentityMatrix[3];
consMat];
ConstraintMatrixC[body_,Hinge]:=Module[{inbAngVelColumns,inbVelColumns,bodyAngVelColumns,bodyVelColumns,inbTrans,consBasis,consMat},
(* cols in constraint matrix corresponding to the inboard body: *)
inbAngVelColumns=(* where ground is not a body *)Flatten[If[Inboard[body]===ground,{{},{},{}}, (* else *)6(Position[Bodies,Inboard[body]][[1]][[1]]-1)]+{1,2,3}];
inbVelColumns=inbAngVelColumns+3; (* translational velocities are right after ang vels *)
bodyAngVelColumns=6(Position[Bodies,body][[1]][[1]]-1)+{1,2,3};
bodyVelColumns=bodyAngVelColumns+3;
inbTrans = CastMtx[Inboard[body],body];(* transform from inboard to body *)
consMat=ConstantArray[0,{5,6Length[Bodies]}]; (* 5 x N array *)
(* Constrain hinge joint translation: *)
consMat[[{1,2,3},inbVelColumns]]=inbTrans; (* velocity of inboard COM *)
       consMat[[{1,2,3},inbAngVelColumns]]=inbTrans.Crossify[BodyTransform[-InbToJnt[body],Inboard[body]]]; (* omega cross r, for inboard COM to joint *)
        consMat[[{1,2,3},bodyVelColumns]]=-IdentityMatrix[3];(* -velocity of body COM *)consMat[[{1,2,3},bodyAngVelColumns]]=-Crossify[BodyTransform[-BodyToJnt[body],body]]; (* - omega cross r, for body COM to joint *)
(* Next constrain hinge joint rotation *)
consBasis=NullSpace[{BodyTransform[JntAxes[body],body]}];(* 2 vectors normal to joint axis *)
consMat[[{4,5},inbAngVelColumns]]=consBasis.inbTrans; (* (omega1 - omega2) dot axis = 0 *)
consMat[[{4,5},bodyAngVelColumns]]=-consBasis;
consMat];
ConstraintMatrixC[body_,Slider]:=Module[{inbAngVelColumns,inbVelColumns,bodyAngVelColumns,bodyVelColumns,inbTrans,consMat,consBasis},
(* cols in constraint matrix corresponding to the inboard body: *)inbAngVelColumns=(* where ground is not a body *)Flatten[If[Inboard[body]===ground,{{},{},{}}, (* else *)6(Position[Bodies,Inboard[body]][[1]][[1]]-1)]+{1,2,3}];
inbVelColumns=inbAngVelColumns+3; (* translational velocities are right after ang vels *)
bodyAngVelColumns=6(Position[Bodies,body][[1]][[1]]-1)+{1,2,3};
bodyVelColumns=bodyAngVelColumns+3;
inbTrans = CastMtx[Inboard[body],body];(* transform from inboard to body *)
consMat=ConstantArray[0,{5,6Length[Bodies]}]; (* 5 x N array *)
(* Constrain slider joint translation *)
       consBasis=NullSpace[{BodyTransform[JntAxes[body],body]}];
(* 2 vectors normal to slider joint axis *)
(* CB (vci + wi x itj = vcb + wb x (btj - jtj) - u axis), where u axis is in nullspace of constraints and thus disappears. *)
consMat[[{1,2},inbVelColumns]]=consBasis.inbTrans; (* velocity of inboard COM *)
       consMat[[{1,2},inbAngVelColumns]]=consBasis.inbTrans.Crossify[BodyTransform[-InbToJnt[body],Inboard[body]]]; 
(* ^ omega cross r, for inboard COM to joint *)
        consMat[[{1,2},bodyVelColumns]]=-consBasis;(* -velocity of body COM *)consMat[[{1,2},bodyAngVelColumns]]=-consBasis.Crossify[BodyTransform[-BodyToJnt[body]+JntToJnt[body],body]]; (* - omega cross r, for body COM to joint *)
(* Next constrain slider joint rotation (none allowed) *)
consMat[[{3,4,5},inbAngVelColumns]]=-inbTrans; (* (omega1 - omega2) = 0 *)
consMat[[{3,4,5},bodyAngVelColumns]]=IdentityMatrix[3];
consMat];
ConstraintMatrixC[body_,UJoint]:=Module[{inbAngVelColumns,inbVelColumns,bodyAngVelColumns,bodyVelColumns,inbTrans,consBasis,consMat},
(* cols in constraint matrix corresponding to the inboard body: *)
inbAngVelColumns=(* where ground is not a body *)Flatten[If[Inboard[body]===ground,{{},{},{}}, (* else *)6(Position[Bodies,Inboard[body]][[1]][[1]]-1)]+{1,2,3}];
inbVelColumns=inbAngVelColumns+3; (* translational velocities are right after ang vels *)
bodyAngVelColumns=6(Position[Bodies,body][[1]][[1]]-1)+{1,2,3};
bodyVelColumns=bodyAngVelColumns+3;
inbTrans = CastMtx[Inboard[body],body];(* transform from inboard to body *)
consMat=ConstantArray[0,{4,6Length[Bodies]}]; (* 4 x N array *)
(* Constrain UJoint translation: *)
consMat[[{1,2,3},inbVelColumns]]=inbTrans; (* velocity of inboard COM *)
       consMat[[{1,2,3},inbAngVelColumns]]=inbTrans.Crossify[BodyTransform[-InbToJnt[body],Inboard[body]]]; (* omega cross r, for inboard COM to joint *)
        consMat[[{1,2,3},bodyVelColumns]]=-IdentityMatrix[3];(* -velocity of body COM *)consMat[[{1,2,3},bodyAngVelColumns]]=-Crossify[BodyTransform[-BodyToJnt[body],body]]; (* - omega cross r, for body COM to joint *)
(* Next constrain UJoint rotation *)
consBasis=NullSpace[{BodyTransform[JntAxes[body],body]}];(* 2 vectors normal to joint axis *)
consMat[[{4},inbAngVelColumns]]=consBasis.inbTrans; (* (omega1 - omega2) dot axis = 0 *)
consMat[[{4},bodyAngVelColumns]]=-consBasis;
consMat];
ConstraintMatrixC[body_,Ball|Gimbal]:=Module[{inbAngVelColumns,inbVelColumns,bodyAngVelColumns,bodyVelColumns,inbTrans,consBasis,consMat,posn},
(* cols in constraint matrix corresponding to the inboard body: *)
inbAngVelColumns=(* where ground is not a body *)Flatten[If[Inboard[body]===ground,{{},{},{}}, (* else *)6(Position[Bodies,Inboard[body]][[1]][[1]]-1)]+{1,2,3}];
inbVelColumns=inbAngVelColumns+3; (* translational velocities are right after ang vels *)
bodyAngVelColumns=6(Position[Bodies,body][[1]][[1]]-1)+{1,2,3};
bodyVelColumns=bodyAngVelColumns+3;
inbTrans = CastMtx[Inboard[body],body];(* transform from inboard to body *)
consMat=ConstantArray[0,{3,6Length[Bodies]}]; (* 3 x N array *)
(* Constrain hinge joint translation *)
       consMat[[{1,2,3},inbVelColumns]]=inbTrans; (* velocity of inboard COM *)
       consMat[[{1,2,3},inbAngVelColumns]]=inbTrans.Crossify[BodyTransform[-InbToJnt[body],Inboard[body]]]; (* omega cross r, for inboard COM to joint *)
        consMat[[{1,2,3},bodyVelColumns]]=-IdentityMatrix[3];(* -velocity of body COM *)consMat[[{1,2,3},bodyAngVelColumns]]=-Crossify[BodyTransform[-BodyToJnt[body],body]]; (* - omega cross r, for body COM to joint *)
(* Ball joint has no rotation constrants, so we're done *)
consMat];


(* ::Section:: *)
(*Generalized Coordinates and Speeds*)


(* ::Input:: *)
(*\[AliasDelimiter]*)


u[n_] := ut[n][t];
If[ TrueQ[ $VersionNumber >= 3.],
	Format[ ut[n_][t] ] := Subscript[u,n];
	Format[ ut[n_]'[t] ] := Subscript[u,n]',
(*else *)
	Format[ ut[n_][t] ] := HoldForm[ Subscripted[u[n]] ];
	Format[ ut[n_]'[t] ] := HoldForm[ Subscripted[u[n]]']
]
	
SetAttributes[ u, Listable ];
attr = If[ TrueQ[$VersionNumber >= 3.],
	NHoldAll, NProtectedAll];
SetAttributes[ u, attr];



q[n_] := qt[n][t];
If[ TrueQ[ $VersionNumber >= 3.],
	Format[ qt[n_][t] ] := Subscript[q,n];
	Format[ qt[n_]'[t] ] := Subscript[q,n]';
	Format[ qt[n_]''[t] ] := Subscript[q,n]'',
(* else *)
	Format[ qt[n_][t] ] := HoldForm[ Subscripted[q[n]] ];
	Format[ qt[n_]'[t] ] := HoldForm[ Subscripted[q[n]]'];
	Format[ qt[n_]''[t] ] := HoldForm[ Subscripted[q[n]]'']
];
SetAttributes[ q, Listable ];
SetAttributes[ q, attr ];


u[n_]' := ut[n]'[t];
q[n_]' := qt[n]'[t];
q[n_]'' := qt[n]''[t];


If[ TrueQ[ $VersionNumber >= 3.],
	Format[ ut[n_][t_] ] := Subscript[u,n][t];
	Format[ qt[n_][t_] ] := Subscript[q,n][t];
	Format[ ut[n_]'[t_] ] := Subscript[u,n]'[t];
	Format[ qt[n_]'[t_] ] := Subscript[q,n]'[t];
	Format[ qt[n_]''[t_] ] := Subscript[q,n]''[t],
(* else *)
	Format[ ut[n_][t_] ] := HoldForm[Subscripted[u[n]][t]];
	Format[ qt[n_][t_] ] := HoldForm[Subscripted[q[n]][t]];
	Format[ ut[n_]'[t_] ] := HoldForm[Subscripted[u[n]]'[t]];
	Format[ qt[n_]'[t_] ] := HoldForm[Subscripted[q[n]]'[t]];
	Format[ qt[n_]''[t_] ] := HoldForm[Subscripted[q[n]]''[t]]
];


u[n_][t_] := ut[n][t]
q[n_][t_] := qt[n][t]
u[n_]'[t_] := ut[n]'[t]
q[n_]'[t_] := qt[n]'[t]
q[n_]''[t_] := qt[n]''[t]


SetAttributes[ ut, Listable ]
SetAttributes[ qt, Listable ]
SetAttributes[ ut, attr ];
SetAttributes[ qt, attr ];


udots[udofs_List] := Map[ ut[#]'[t]&, udofs ]
udots[udof_] := ut[udof]'[t]
qdots[qdofs_List] := Map[ qt[#]'[t]&, qdofs ]
qdots[qdof_] := qt[qdof]'[t]


States := Join[ q[qdofs], u[udofs] ];


(* ::Section:: *)
(*Positions on Bodies*)


PosCOM[ body_ /; MemberQ[Bodies, body] ] :=
  PosCOM[ Inboard[body] ] + InbToJnt[body]+JntToJnt[body]-BodyToJnt[body];


(* ::Input:: *)
(*PosCOM[ body_  ] :=*)
(*  PosCOM[ Inboard[body] ] + InbToJnt[body]+JntToJnt[body]-BodyToJnt[body];*)


PosCOM[(bodies_)?ListQ] := 
  Simplify[Plus @@ (Mass /@ bodies*PosCOM /@ bodies)/
    Plus @@ Mass /@ bodies]


PosPnt[ point_, body_ ] := PosCOM[ body ] + point;


Dist[ pnt1_, pnt2_ ] :=
	Simplify[ NormV[ Transpose[ Level[
		CastV[ pnt2 - pnt1 , ground ], 1]][[1]] ] ]


(* ::Section:: *)
(*Velocities on Bodes*)


VelPntFixed[point_, body_] := VelCOM[body] + Cross[AngVel[body], point]; 
AccPntFixed[point_, body_] := AccCOM[body] + Cross[AngAcc[body], point] + Cross[Cross[AngVel[body], point]]; 


(* ::Section:: *)
(*Derivative of Free Vector*)


(* ::Special1:: *)
(*Partial and total derivatives are computed separately. Partial derivatives of single or multiple vectors are found using D[ vector, scalar, in->frame, NonConstants->{functions} ], where in and NonConstants are options. Option in is used to specify the reference frame the partial derivative is desired in, and defaults to ground.*)
(**)
(*Total derivatives of a single or multiple vector are found using Dt[ vector, scalar, in->frame, Constants->{constants} ], where in and Constants are options. Option in is used to specify the reference frame the total derivative is desired in, and defaults to ground. If the angular velocities of all of the reference frames vector is defined in are themselves defined, the total time derivative is computed using those angular velocities.*)


Options[DPV] = { in->ground };


PV/: D[ PV[{c_,frme_,n_}], x_, opts___Rule ] := Module[
	{bfrme, ncons},
	bfrme = in  /. {opts} /. Options[DPV];
	ncons = NonConstants /. {opts} /. NonConstants->{};
	Dprtl[ PV[{c,frme,n}], x, bfrme, NonConstants->ncons ] ];
PV/: D[ PV[vecs:{_,_,_}..], x_, opts___Rule] := Apply[ Plus,
	Map[ D[ PV[#], x, opts]&, {vecs} ] ];


PV/: Dprtl[ PV[{c_,frme_,n_}], t, bfrme_, opts___Rule] :=
	Module[ {pv},
		pv = Apply[ List, CastV[ PV[{c,frme,n}], bfrme ] ] /.
			{q[k_] -> qh[k], u[k_] -> uh[k]};
		Apply[ PV, Transpose[
			MapAt[ D[#,x,opts]&, Transpose[pv], 1] ] ] /.
				{qh[k_]->q[k],uh[k_]->u[k]} ];


PV/: Dprtl[ PV[{c_,frme_,n_}], x_, bfrme_, opts___Rule] :=
	Module[ {pv},
		pv = Apply[ List, CastV[ PV[{c,frme,n}], bfrme ] ];
		Apply[ PV, Transpose[
			MapAt[ D[#,x,opts]&, Transpose[pv], 1] ] ] ];


PV/: Dt[ PV[{c_,frme_,n_}], x_, opts___Rule] := Module[
	 {bfrme, cons},
	 bfrme = in /. {opts} /. Options[DPV];
	 cons = Constants /. {opts} /. Constants->{};
	 Dttl[ PV[{c,frme,n}], x, bfrme, Constants->cons ] ];


PV/: Dt[ PV[v:{_,_,_}..], x_, opts___Rule] :=
	  Apply[ Plus, Map[ Dt[ PV[#], x, opts]&, {v} ] ];


PV/: Dttl[ PV[{c_,frme_,n_}], t, bfrme_, opts___Rule] :=
	Cross[ AngVel[frme,bfrme], PV[{c,frme,n}] ] +
	PV[{Dt[c,t,opts],frme,n}] /; 
		ValueQ[AngVel[frme]] && ValueQ[AngVel[bfrme]]


PV/: Dttl[ PV[{c_,frme_,n_}], x_, bfrme_, opts___Rule] :=
	Module[ {pv},
		pv = Apply[ List, CastV[ PV[{c,frme,n}], bfrme ] ];
		Apply[ PV, Transpose[
			MapAt[ Dt[#,x,opts]&, Transpose[pv], 1] ] ] ];


(* ::Input:: *)
(*(* Older versions of derivatives *)*)


(* ::Input:: *)
(*Null*)


(* ::Input:: *)
(*PV/: Dt[ PV[{c_,frme_,n_}], t, opts___] := *)
(*  Cross[ AngVel[frme], PV[{c,frme,n}] ]+*)
(*  PV[{Dt[c,t,opts],frme,n}];*)
(*PV/: Dt[ PV[x:{_,_,_}..], t, opts___] :=*)
(*  Apply[ Plus, Map[ Dt[ PV[#], t, opts]&, {x} ] ];*)


(* ::Input:: *)
(*PV/: D[ PV[{c_,frme_,n_}], t] :=*)
(*  Cross[ AngVel[frme], PV[{c,frme,n}] ]+*)
(*  PV[{D[c,t],frme,n}];*)


(* ::Input:: *)
(*PV/: D[ PV[x:{_,_,_}..], t] :=*)
(*	Apply[ Plus, Map[ D[ PV[#], t]&, {x} ] ];*)


(* ::Section:: *)
(*Angular Momentum*)


AngMom[ body_, pnt_ ] :=
	Inertia[body].AngVel[body] +
	Cross[ Mass[body] (PosCOM[body] - pnt),
		VelCOM[body] ];
AngMom[ bodies_List, pnt_ ] :=
	Apply[ Plus, Map[ AngMom[#, pnt]&, bodies ] ];
AngMom[ body_, pnt_, ref_ ] :=
	Inertia[body].AngVel[body,ref] +
	Cross[ Mass[body] (PosCOM[body] - pnt),
		VelCOM[body,ref] ];
AngMom[ bodies_List, pnt_, ref_ ] :=
	Apply[ Plus, Map[ AngMom[#, pnt, ref]&, bodies ] ];
AngMom[ body_ ] := Inertia[body].AngVel[body]; (* body about its own com *)


(* ::Section:: *)
(*Linear Momentum*)


LinMom[body_] := Mass[body]*VelCOM[body]; 
  LinMom[(bodies_)?ListQ] := Simplify[Plus @@ LinMom /@ bodies]; 


(* ::Section:: *)
(*Partial velocities*)


(* ::Special1:: *)
(*To compute partial velocities, you need to be able to take the partial derivative of an expression with respect to a u.*)


PV/: D[ PV[ {c_, frme_, n_ }], u[d_] ] :=
  PV[{ D[c,u[d]], frme, n}];
PV/: D[ PV[x:{_,_,_}..], u[d_] ] := 
	Apply[Plus, Map[ D[PV[#], u[d]]&, {x}] ];


PrtVel[ expr_, n_?NumberQ ] := D[ expr, u[n] ];
PrtVel[ expr_, u[n_?NumberQ] ] := D[ expr, u[n] ];


(* ::Section:: *)
(*Generalized Active Forces*)


ResltF[ body_ ] := Apply[ Plus, Map[ #[[1]]&, Force[body] ] ];
ResltT[ body_ ] := Apply[ Plus, Torque[body] ] + 
	Apply[ Plus, Map[ Cross[ #[[2]], #[[1]] ] &, Force[body] ] ];
GActFrcC[ All ]  := (* Simplify[ *)
	Apply[ Plus, Map[ GActFrc, Bodies ] ] (*]*);
GActFrcC[ All, n_ ] := Apply[ Plus,
	Map[ GActFrc[#,n]&, Bodies] ];
GActFrcC[ body_, udof_ /; MemberQ[udofs,udof] ] :=
		PrtVel[ AngVel[body], udof] . ResltT[body] +
		PrtVel[ VelCOM[body], udof] . ResltF[body] /.
		Join[replacesin,replacecos,replace2SC,replaceSC2] //ZF;
GActFrcC[ body_ ] := Map[ GActFrc[ body, # ]&, udofs ];
SetAttributes[ GActFrcC, Listable ]


AppFrc[ body_, force_, point_ ] := (
	AppendTo[ Force[body], {force, point} ];
	Length[ Force[body] ] )
AppTrq[ body_, torque_ ] := (
	AppendTo[ Torque[body], torque ];
	Length[ Torque[body] ] )


(* ::Section:: *)
(*Generalized Inertia Forces*)


(* ::Input:: *)
(*(* Note: C stands for conventional, S stands for simplified, 0 stands for bare minimum *)*)


GInerFrcC[ All ] := CollectU[ 
	Apply[ Plus, Outer[ GInerFrc, Bodies, udofs ] ] ] //.
		Join[replacesin,replacecos,replace2SC,replaceSC2];
	(* Apply[ Plus, Map[ gInerFrc, Bodies ] ]; *)
GInerFrcC[ All, n_ ] := CollectU[
	Apply[ Plus, Map[ GInerFrc[#,n]&, Bodies ] ] ] //.
		Join[replacesin,replacecos,replace2SC,replaceSC2];
GInerFrcC[ body_, udof_ /; MemberQ[udofs,udof] ] :=
	Collect[
		PrtVel[ AngVel[body], udof] . Tstar[body] +
		PrtVel[ VelCOM[body], udof] . Rstar[body] //.
		Join[replacesin,replacecos,replace2SC,replaceSC2],
		Join[ Map[(u[#])'&, udofs], Map[u[#]&, udofs] ] ] //ZF
GInerFrcC[ body_ ] := CollectU[
	Map[ GInerFrc[ body, # ]&, udofs ] ] //.
		Join[replacesin,replacecos,replace2SC,replaceSC2];
GInerFrc0[ All ] := 
	Apply[ Plus, Outer[ GInerFrc, Bodies, udofs ] ]  //.
		Join[replacesin,replacecos,replace2SC,replaceSC2];
	(* Apply[ Plus, Map[ GInerFrc, Bodies ] ]; *)
GInerFrc0[ All, n_ ] := 
	Apply[ Plus, Map[ GInerFrc[#,n]&, Bodies ] ]  //.
		Join[replacesin,replacecos,replace2SC,replaceSC2];
GInerFrc0[ body_, udof_ /; MemberQ[udofs,udof] ] :=
	Collect[
		PrtVel[ AngVel[body], udof] . Tstar[body] +
		PrtVel[ VelCOM[body], udof] . Rstar[body] //.
		Join[replacesin,replacecos,replace2SC,replaceSC2],
		Join[ Map[(u[#])'&, udofs], Map[u[#]&, udofs] ] ] //ZF
GInerFrc0[ body_ ] := 
	Map[ GInerFrc[ body, # ]&, udofs ] //.
		Join[replacesin,replacecos,replace2SC,replaceSC2];

GInerFrcS[ All ] := CollectU[ 
	Apply[ Plus, Outer[ GInerFrc, Bodies, udofs ] ] ] //.
		Join[replacesin,replacecos,replace2SC,replaceSC2];
GInerFrcS[ All, n_ ] := CollectU[
	Apply[ Plus, Map[ GInerFrc[#,n]&, Bodies ] ] ] //.
		Join[replacesin,replacecos,replace2SC,replaceSC2];
GInerFrcS[ body_ ] := CollectU[
	Map[ GInerFrc[ body, # ]&, udofs ] ] //.
		Join[replacesin,replacecos,replace2SC,replaceSC2];
GInerFrcS[ body_, udof_ /; MemberQ[udofs,udof] ] :=
	Simplify[Collect[
		PrtVel[ AngVel[body], udof] . Tstar[body] +
		PrtVel[ VelCOM[body], udof] . Rstar[body] //.
		Join[replacesin,replacecos,replace2SC,replaceSC2],
		Join[ Map[(u[#])'&, udofs], Map[u[#]&, udofs] ] ] ] //ZF

SetAttributes[ GInerFrcC, Listable ]
SetAttributes[ GInerFrcS, Listable ]
SetAttributes[ CleanGIF, Listable ]
SetAttributes[ CleanGIFS, Listable ]


TstarC[body_] := (-AngAcc[body].Inertia[body] - 
	Cross[ AngVel[body], Inertia[body] . AngVel[body] ]) //ZF
Rstar[body_] := -Mass[body] AccCOM[body]
SetAttributes[ TstarC, Listable ]
SetAttributes[ Rstar, Listable ]


CollectU[ x_ ] := Collect[ x, 
	Join[ Map[(u[#])'&, udofs], Map[u[#]&, udofs] ] ];


upow2 := Union[ Flatten[ 
	Outer[ Times, 
		Map[u[#]&, udofs], Map[u[#]&, udofs] 
	] ] ]


CleanGIF[ x_ ] := Apply[ Plus, Map[ 
	Coefficient[Expand[x], #] # &, 
	Join[ Map[(u[#])'&, udofs], upow2 ] ] ]


CleanGIFS[ x_ ] := Apply[ Plus, Map[ 
	Simplify[ Coefficient[Expand[x], #] ] # &, 
	Join[ Map[(u[#])'&, udofs], upow2 ] ] ]


(* ::Section:: *)
(*Zees*)


SetAttributes[Zee,{Orderless,Listable}];


Z[n_]:=Zz[n][t];


Format[Zz[n_][t_]]:=Subscript[Z,n];
Format[ZC[n_]]:=Subscript[ZC,n];


(* ::Input:: *)
(*Zee[x_List]:=Zee/@x;*)


(* Find constants independent of t, and make them zee constants *)
Zee[x_?(FreeQ[#,t]&)  y_.]:=Module[{xs=x/.Join[replacesin,replacecos,replace2SC,replaceSC2]},
If[LeafCount[xs]>4,
 ZC[addtozeeC[xs]]y, (* else *)
xs y]];


(* Make coefficients of u' into zees, if above a minimum complexity *)
Zee[x_ ut[n1_]'[t_]]:=Module[{xs=x/.Join[replacesin,replacecos,replace2SC,replaceSC2]},
If[LeafCount[xs]>4,Z[addtozee[xs]] ut[n1]'[t],
		xs ut[n1]'[t] ]];


(* Make coefficients of u*u into zees, if above a minimum complexity *)
Zee[x_ ut[n1_][t_] ut[n2_][t_]]:=Module[{xs=x/.Join[replacesin,replacecos,replace2SC,replaceSC2]},
If[LeafCount[xs]>4,Z[addtozee[xs]] ut[n1][t]ut[n2][t],
		xs ut[n1][t] ut[n2][t]]];


Zee[x_ ut[n1_][t_] ut[n1_][t_]]:=Module[{xs=x/.Join[replacesin,replacecos,replace2SC,replaceSC2]},
If[LeafCount[xs]>4,Z[addtozee[xs]] ut[n1][t]ut[n1][t],
		xs ut[n1][t] ut[n1][t]]];


(* Make coefficients of u into zees, if above a minimum complexity *)
Zee[x_ ut[n1_][t_]  /;FreeQ[x,u[_Integer]]]:=Module[{xs=x/.Join[replacesin,replacecos,replace2SC,replaceSC2]},
If[LeafCount[xs]>4,Z[addtozee[xs]] ut[n1][t],
		xs ut[n1][t] ]];


(* Within a vector, collect things linear in u or u' and make them into zees *)
Zee[PV[List[x_,y_,z_]]]:=PV[{Zee[Collect[x,Join[u[udofs],udots[udofs]]]],y,z}]


Zee[PV[x_]]:=Transpose[Map[Zee[Collect[#,Join[u[udofs],udots[udofs]]]]&,Transpose[x][[1]]],Drop[Transpose[x],1]];


(* For a sum of terms, apply zees to linear coefficients of u and u' *)Zee[x:Plus[__]]:=Map[Zee[Collect[#,Join[u[udofs],udots[udofs]]]]&,x];


(* Be able to handle time derivatives of zees *)
Zz[n_]'[t]:=D[zees[[n]],t];


addtozee[x_ /;!MemberQ[zees,x]]:=Module[{},AppendTo[zees,x];Length[zees]]
addtozee[x_/;MemberQ[zees,x]]:=Module[{},Flatten[Position[zees,x,1]][[1]]]


addtozeeC[x_/; !MemberQ[zeesC,x]]:=Module[{},AppendTo[zeesC,x];Length[zeesC]]
addtozeeC[x_/; MemberQ[zeesC,x]]:=Module[{},Flatten[Position[zeesC,x,1]][[1]]]


UnZee[x_]:=x//.{Z[n_]:>zees[[n]],ZC[n2_]:>zeesC[[n2]]};


(* Given an expression or list of expressions, find exhaustively all the ZCs and Zs referred to within, and follow the definitions for all of these ZCs and Zs to yield a minimum set of Zee indices. Output is in the form of {zeesC, zees} that are used *)
MinimumZee[expr_]:=Module[ (* start with a list of the ZCs and Zs referred to explicitly by the expressions *){startingzees={Union[Cases[expr,ZC[n_]:>n,Infinity]],Union[Cases[expr,Zz[n_][t]:>n,Infinity]]}},
FixedPoint[zeesbothmoreused,startingzees]]


(*Helper function for MinimumZee, is called repeatedly to look at the definitions of the ZCs and Zs and find what other Zs are referred to therein. Attach those to the existing list, and return a new list of ZCs and Zs. *)
zeesbothmoreused[{usedzeecs_List,usedzees_List}]:={
Union[(* identify ZCs referenced by existing ZCs and Zs*)Cases[zeesC[[#]]& /@usedzeecs (* a list of defs for zeesC referred to by existing zeesC *),
ZC[n_]:>n,Infinity] (* Cases extracts the ZC numbers *) ~Join~
Cases[zees[[#]]&/@usedzees (* a list of defs for the zees referred to by existing zees *),
ZC[n_]:>n,Infinity](* Cases extracts the ZC numbers therein *)
~Join~usedzeecs] (* and grow onto the list of ZCs *),
Union[ (* identify zees referenced by existing Zs *)Cases[zees[[#]]&/@usedzees(* a list of defs for zees referred to by existing zees *),
Zz[n_][t]:>n,Infinity](* Cases extracts the Z numbers therein *)
~Join~usedzees]};



(* ::Section:: *)
(*Equations of Motion*)


Options[EOM] = {Simplify -> Off}; 


EOM[opts___] := Module[{eom},
eom = MapAll[ZF,
If[(Simplify/.{opts}/.Options[EOM])===On,
  (* then Simplify it *)
  Thread[ 
	GActFrc[All] + GInerFrcS[All] == 0*udofs ],
  (* else regular equations *)
  Thread[ 
	GActFrc[All] + GInerFrc[All] == 0*udofs ]
] ];
(* One last attempt to identify constant factors in zees *)
zees = ZF[zees]; (* This should only affect zeesC *)
eom ]

(* EOM0 := Thread[ 
	GActFrc[All] + GInerFrc0[All] == Table[0, {Length[udofs]} ] ] *)


(* ::Special1:: *)
(*CollectE takes the equations of motion, and collects the terms of the udots together. You can also use CollectE to collect terms linear in a specified variable x, with CollectE[ eom, x].*)


CollectE[ eom_ ] := Collect[eom,
	udots[udofs] ] == 0  
CollectE[ eom_, x_ ] := Collect[eom,
	x ] == 0
SetAttributes[ CollectE, Listable ]


(* ::Special1:: *)
(*SolveEqn solves a list of equations for the given variables, sorts them into order, and puts an == sign around the whole thing*)


SolveEqn[ eqn_, vars_ ] := Sort[ Apply[ Equal, 
	Solve[ eqn, vars], {2}][[1]] ]


(* ::Special1:: *)
(*RHS gives the right-hand-side of each equation in a list of equations*)


RHS[ eqn_ ] := Map[#[[2]]&,eqn]


(* ::Special1:: *)
(*AsmStateEqn assembles the state equations, if you give it the kinematical differential equations and the equations of motion. First, it solves them for the q's and u's, then it joins it all together, and provides the right-hand-side. This output is suitable for RungeKutta*)


AsmStateEqn[ Kinematics_, eom_ ] := RHS[ Join[ 
	SolveEqn[ Kinematics, qdots[qdofs] ],
	SolveEqn[ eom, udots[udofs] ] ]]


(* ::Special1:: *)
(*SolveKinematics gives a solution to the kinematical equations, to allow qdots to be replaced by u' s.*)


SolveKinematics:=Solve[Kinematics,qdots[qdofs]]//Flatten


(* ::Special1:: *)
(*MassMatrix extracts the mass matrix from the equations of motion*)


MassMatrix[ eom_ ] := Module[ {eqns,mass},
	eqns = Thread[ eom, Equal]; eqns = eqns[[1]]-eqns[[2]];
	mass = Transpose[Map[ Coefficient[ Expand[eqns], ut[#]'[t] ] &, udofs ]];
  {-mass, (eqns /. Thread[udots[udofs]->Table[0,{Length[udofs]}]])}
] 


(* ::Special1:: *)
(*XDot computes the state derivative, given the equations of motion and the state*)


XDot[eom_, states_List] := Module[
	{staterule = Thread[
		Evaluate[ Join[q[qdofs],u[udofs]] ] -> 
			Flatten[states] ],
	xdots = 
		Join[Thread[q[qdofs]'],Thread[u[udofs]'] ]},
	xdots /. Solve[
		Flatten[{eom /. parms, Kinematics}] /.
		staterule, xdots][[1]]
]


(* ::Section:: *)
(*Exporting*)


lowercasefunctions = {Sin -> sin, Cos -> cos, Tan -> tan, Power -> power,Sqrt->sqrt,
ArcSin->asin,ArcCos->acos,ArcTan[x_]->atan[x],ArcTan[x_,y_]->atan2[y,x]}; 


(* ::Text:: *)
(*formtrigreplacements takes in a set of equations, and returns a list of replacements*)


Options[formtrigreplacements]={TrigReplacements->{{Sin,"s"},{Cos,"c"},{Tan,"t"},{Sec,"sec"},{Csc,"csc"}},ArgReplacements->None };
(* trigSubs = formtrigreplacements[equations] takes in equations of any form, and returns a list of replacements for trig function calls with shorthand variables. The shortened variables are given as strings. For example, an equation M*g*Sin[q[1]]*Cos[q[2]]*Sin[q[1]-q[2]]*Sin[unknown] returns {Sin[q[1]]->"s1",Cos[q[2]]->"c2",Sin[q[1]-q[2]]->"s1m2"}. There are two optional arguments that can be given. ArgReplacements\[Rule]{{gamma,"g"},{alpha,"a"}} specifies additional trig arguments that can be shortened, and TrigReplacements\[Rule]{{Sin,"s"},{Cos,"c"}} specifies all the trig functions that can be shortened. The second element of each entry is a string containing the abbreviation.  *)
formtrigreplacements[eom_,opts___]:=Module[{trigreps,altargs,alltrig,trigreplaced,argreplacements,argtostrings,plusminustostrings,successfulones,trigreplacements,trigfncalls,trigSubs,reorder},
trigreps=TrigReplacements/.{opts}/.Options[formtrigreplacements]/.None->{};
altargs=ArgReplacements/.{opts}/.Options[formtrigreplacements]/.None->{};
	(* pull out all the unique trig function calls of the form Sin[_]|Cos[_] etc., yielding a list like {Sin[q[1]],Cos[q[2]-q[3]]} *)
alltrig=Union[Cases[eom,Apply[Alternatives,Map[#[_]&,Transpose[trigreps][[1]]]],Infinity]];
	(* form a list of the trig functions replaced by their shortcuts, e.g. {Cos[q[1]],Sin[anything]} becomes {"c","s"} *)
trigreplaced=Map[Head,alltrig]/.Apply[Rule,trigreps,{1}];
	(* make a list of replacements, starting with q[n_]\[RuleDelayed]ToString[n], followed by alternative arg replacements if any (whether in vector or matrix form) e.g. {q[n_]\[RuleDelayed]ToString[n],gamma->"g"} *)
argreplacements=Prepend[Apply[Rule,Replace[altargs/.None->{},x:{_,_String}->{x}],{1}],q[n_]:>ToString[n]];
	(* convert all arguments into their replacements, e.g. q[1] becomes "1" *)
argstostrings=Map[#[[1]]&,alltrig]/.argreplacements;
	(* for cases where two arguments are summed, convert all plus and minus signs into "p" and "m" *)
plusminustostrings=argstostrings/.(x_String - y_String):>(x <> "m" <> y)/.(x_String + y_String) :>( x <> "p" <> y);
	(* positions of the arguments that were successfully replaced *)
successfulones=Position[plusminustostrings,_String,{1}]; 
	(* splice together the trig function replacement (e.g. "s" for Sin) with the argument replacements, yielding something like {"s1","c2","s12"} *)
trigreplacements=MapThread[StringJoin,{Extract[trigreplaced,successfulones],Extract[plusminustostrings,successfulones]}];
(* and pull out the matching trig functions, effectively discarding the ones that couldn't be replaced, e.g. {Sin[q[1],Cos[q[2]],Sin[q[1]-q[2]]} *)
(* sort them so that shorter strings go first, and after that they are alphabetical *)
reorder=Ordering[trigreplacements,All,StringLength[#1]<StringLength[#2]||(StringLength[#1]===StringLength[#2] && OrderedQ[{#1,#2}])&];
trigfncalls=Extract[alltrig,successfulones];
(* Use the two lists to form a set of rules, {Sin[q[1]]\[Rule]"s1",Cos[q[2]]->"c2"} *)
trigSubs=MapThread[Rule[#1,#2]&,{trigfncalls[[reorder]],trigreplacements[[reorder]]}]]


formqureplacements:=Module[{qSubs,uSubs,upSubs,quSubs,activeqdofs,activeudofs},
activeqdofs=Select[qdofs,q[#]=!= 0&]; (* deactive qs and us by setting them to zero *)
activeudofs=Select[udofs,u[#]=!=0&];
qSubs=(q[#1]->ToString[StringForm["q`1`", #]]&)/@activeqdofs;uSubs=(u[#1]->ToString[StringForm["u`1`", #]]&)/@activeudofs;upSubs=(Derivative[1][u[#1]]->ToString[StringForm["u`1`dot", #]]&)/@activeudofs;
quSubs=Join[qSubs,uSubs,upSubs]];


(* ::Text:: *)
(*We will change anything raised to the 2nd power into a multiplication, because Mathematica wants to export squares with the Power function*)


pSub = {(x_)^2 -> HoldForm[x*x]}; 


(* ::Text:: *)
(*ExportExpression[expression, "name"] formats expression as Matlab code, producing a string that begins with 'name = ' followed by the expression. ExportExpression can export a matrix, a symmetric matrix, a list, or a scalar. Formatting options include EndLine->"\n", Continuator->"...", LineLength->78, Deliminator->" ", which specify the end-of-line character, line continuator, maximum line length of the output, and the deliminators where line breaks can occur (default to space character). The options to formtrigreplacements can also be used. Also: ExportExpression[expr, "name", Indices->indices] exports only the listed indices within expr.*)


Options[ExportExpression]={qusubs->Automatic,TrigSubs->Automatic,psub->{(x_)^2->HoldForm[x*x]},zsub->{Zz[n_][t]->HoldForm[Z[n]]},Indices->All,Preallocate->True};
(* if a symmetric matrix, export with that knowledge *)
ExportExpression[Mx_?SymmetricMatrixQ, name_String,opts___]:=Module[{elementsreplaced,quSubs,trigSubs,pSub,zSub,preallocation={}},
quSubs=qusubs/.{opts}/.Options[ExportExpression]/.Automatic->formqureplacements;
trigSubs=TrigSubs/.{opts}/.Options[ExportExpression]/.None->{}/.Automatic->formtrigreplacements[Mx,opts];
pSub=psub/.{opts}/.Options[ExportExpression];
zSub=zsub/.{opts}/.Options[ExportExpression]; 
(* produce table of values, such as {{MM(1,1)=M*g*s12,MM(1,2)=c12},{MM(2,1)=c12,MM(2,2)=c2}} for
symmetric matrices, where duplicate definitions are avoided *)
If[Preallocate/.{opts}/.Options[ExportExpression],(* true if want to preallocate *)preallocation=ToString[StringForm["`1` = zeros(`2`,`3`);\n",name,Sequence@@Dimensions[Mx]]]];
elementsreplaced=Table[
If[Mx[[i,j]]===0,"",
StringReplace[ (* strip quotes *)
(* expression of the form name(i,j) = rhs *)
ToString[StringForm["`1`(`2`,`3`) = `4`; ",name,i,j,If[i>j,StringForm["`1`(`2`,`3`)",name,j,i],CForm[Mx[[i,j]]/.trigSubs/.quSubs/.pSub/.zSub]//.lowercasefunctions]]],"\""->""]],
{i,1,Dimensions[Mx][[1]]},{j,1,Dimensions[Mx][[2]]}];
(* output each row as a separate formatted output, with line breaks in between *)
	StringJoin[preallocation,Map[StringInsert[fmtString[#,opts],"\n",-1]&,Map[StringJoin,elementsreplaced]]]
];
(* if a matrix, export accordingly *)
ExportExpression[Mx_List?MatrixQ, name_String,opts___]:=Module[{elementsreplaced,quSubs,trigSubs,pSub,zSub,preallocation={}},
quSubs=qusubs/.{opts}/.Options[ExportExpression]/.Automatic->formqureplacements;
trigSubs=TrigSubs/.{opts}/.Options[ExportExpression]/.None->{}/.Automatic->formtrigreplacements[Mx,opts];
pSub=psub/.{opts}/.Options[ExportExpression];
zSub=zsub/.{opts}/.Options[ExportExpression];
(* produce table of values, such as {{MM(1,1)=M*g*s12,MM(1,2)=c12},{MM(2,1)=c12,MM(2,2)=c2}} *)
If[Preallocate/.{opts}/.Options[ExportExpression],(* true if want to preallocate *)preallocation=ToString[StringForm["`1` = zeros(`2`,`3`);\n",name,Sequence@@Dimensions[Mx]]]];
elementsreplaced=Table[
If[Mx[[i,j]]===0,"",(* else write it out *)
StringReplace[ (* strip quotes *)
(* expression of the form name(i,j) = rhs *)
ToString[StringForm["`1`(`2`,`3`) = `4`; ",name,i,j,CForm[Mx[[i,j]]/.trigSubs/.quSubs/.pSub/.zSub]//.lowercasefunctions]],"\""->""]],
{i,1,Dimensions[Mx][[1]]},{j,1,Dimensions[Mx][[2]]}];(* output each row as a separate formatted output, with line breaks in between *)
	StringJoin[preallocation,Map[StringInsert[fmtString[#,opts],"\n",-1]&,Map[StringJoin,elementsreplaced]]]
];
(* if a vector, export accordingly *)
ExportExpression[rhs_List?VectorQ, name_String,opts___]:=Module[{elementsreplaced,quSubs,trigSubs,pSub, zSub,indices,preallocation={}},
quSubs=qusubs/.{opts}/.Options[ExportExpression]/.Automatic->formqureplacements;
trigSubs=TrigSubs/.{opts}/.Options[ExportExpression]/.None->{}/.Automatic->formtrigreplacements[rhs,opts];
pSub=psub/.{opts}/.Options[ExportExpression];
zSub=zsub/.{opts}/.Options[ExportExpression]; 
indices=Indices/.{opts}/.Options[ExportExpression]/.All->Length[rhs];
(* produce list of values, such as {rhs(1)=M*g*s12,rhs(2)=c12} *)
If[Preallocate/.{opts}/.Options[ExportExpression],(* true if want to preallocate *)preallocation=ToString[StringForm["`1` = zeros(`2`,1);\n",name,Sequence@@Dimensions[rhs]]]];

elementsreplaced=Table[
StringReplace[ (* strip quotes *)
(* expression of the form name(i,j) = rhs *)
ToString[StringForm["`1`(`2`) = `3`; ",name,i,CForm[rhs[[i]]/.trigSubs/.quSubs/.pSub/.zSub]//.lowercasefunctions]],"\""->""],
{i,indices}];
(* output each row as a separate formatted output, with line breaks in between *)
StringJoin[preallocation,Map[StringInsert[fmtString[#,opts],"\n",-1]&,elementsreplaced]]];
(* if a scalar, export accordingly *)
ExportExpression[scalar_ /;Nor[VectorQ[scalar],MatrixQ[scalar] ],name_String ,opts___]:=Module[{quSubs,trigSubs,pSub, zSub},
quSubs=qusubs/.{opts}/.Options[ExportExpression]/.Automatic->formqureplacements;
trigSubs=TrigSubs/.{opts}/.Options[ExportExpression]/.None->{}/.Automatic->formtrigreplacements[scalar,opts];
pSub=psub/.{opts}/.Options[ExportExpression];
zSub=zsub/.{opts}/.Options[ExportExpression]; 
fmtString[StringReplace[ToString[StringForm["`1` = `2`;",name,
CForm[scalar/.trigSubs/.quSubs/.pSub/.zSub]//.lowercasefunctions]],"\""->""],opts]];


(* ::Text:: *)
(*fmtString[string] reformats a long String argument and tries to break it into a bunch of lines with a Matlab (or other) continuator, so that the output is fairly readable code. Options include EndLine->\"\\n\" (end of line character), Continuator->\"...\" (Matlab or other line continuator), LineLength->78 (max line length), Deliminator->\" \" (look for spaces, lists of deliminators are acceptable).*)


Options[fmtString]={EndLine->"\n",Continuator->"...",LineLength->78,Delimiter->" "};
(* fmtString[string] takes a long String argument and tries to break it into a bunch of lines with a Matlab (or other) continuator, so that the output is digestable code. Options include EndLine->"\n" (end of line character), Continuator->"..." (Matlab or other line continuator), LineLength\[Rule]78 (max line length), Delimiter->" " (look for spaces, lists of delimiters are acceptable) *)
fmtString[s_,opts___]:=Module[{maxlen,instring=s,outstring="",spaces,eol,continuator,delim},
eol=EndLine/.{opts}/.Options[fmtString]/.None->"";
continuator=Continuator/.{opts}/.Options[fmtString]/.None->"";
maxlen=LineLength/.{opts}/.Options[fmtString];
delim=Delimiter/.{opts}/.Options[fmtString];
While[StringLength[instring]>maxlen,
spaces=StringPosition[ StringTake[instring, maxlen], delim];
If[ Length[spaces] < 1, spaces = StringPosition[ instring, delim]];
outstring = outstring <> StringTake[instring, Last[spaces][[1]]] <> continuator<>eol;
instring = StringDrop[ instring, Last[spaces][[1]]];
];
outstring = outstring <> instring ]


(* ::Text:: *)
(*LinearTerms[ expressions, terms] returns a matrix of linear coefficients of the variables listed in terms, that appear in the expressions (also a list) given. One option is Simplifier->Simplify, which causes the specified function to be applied to the matrix. UnLinearTerms returns all of the remaining terms that are not linear in the LinearTerms.*)


Options[LinearTerms]={Simplifier->Identity};
LinearTerms[ rhs_,terms_List,opts___ ] := Module[{linearstuff,simp},simp=Simplifier/.{opts}/.Options[LinearTerms];linearstuff=Outer[Coefficient[#1,#2]&,rhs,terms]//simp];
Options[UnLinearTerms]={Simplifier->Identity};
UnLinearTerms[rhs_,terms_,opts___]:=Module[{simp},simp=Simplifier/.{opts}/.Options[UnLinearTerms];
rhs/.Map[#->0&,terms]//simp]
LinearTerms[rhs_,term_,opts___]:=LinearTerms[rhs,{term},opts]; (* if only one term is used, turn it into a list *)


(* ::Text:: *)
(*ExportEOM writes a Matlab function to evaluate the state-derivative. Features include substitution of q's and u's for the state x, and substitution of shorthand trig functions such as s1 for sin(q(1)). ExportEOM[{massmatrix,rhs}] uses the mass matrix and right-hand side given (usually found using MassMatrix[]). Alternatively, ExportEOM[equations] will call MassMatrix[] automatically. Several optional arguments may be given. Forces->{T1[t],T2[t]} specifies forces that enter linearly into the equations; ExportCode will export a matrix of coefficients for these forces. Another option is Expressions->{{KE,"KE"},{PE,"PE"}}, which outputs the specified additional expressions, named with the strings given. Each expression must be paired with a name within a list. The output of ExportEOM can be sent to a file with OutputFile->"filename", which is automatically produced in the current directory. The Matlab function is named after the OutputFile if specified, with default "fstatederivative". It can also be explicitly set with option mFileName->"name". Formatting options include EndLine->"\n", Continuator->"...", LineLength->78, Deliminator->" ", which specify the end-of-line character, line continuator, maximum line length of the output, and the deliminators where line breaks can occur (default to space character). Zees-> All, None, Minimum*)


Options[ExportEOM]={OutputFile->None,Forces->None,Expressions->None,Zees->Minimum,mFileName->"fstatederivative"};


(* ::Input:: *)
(*SetAttributes[ExportEOM2,HoldRest];Options[ExportEOM2]={Expressions->None};*)
(*ExportEOM2[{mass_List?MatrixQ,rhs_List?VectorQ},opts:OptionsPattern[]]:=Module[{unevaled,names,expressions,expressionpairs},*)
(*unevaled=OptionValue[ExportEOM2,Unevaluated[{opts}],Expressions,Hold];*)
(*names=StringSplit[ToString[unevaled],{", ","{"->"","Hold["->"","}]"->"","}"->""}];*)
(*expressions=ReleaseHold[OptionValue[ExportEOM2,{opts},Expressions]];*)
(*expressionpairs=Switch[expressions,*)
(*	None,None,*)
(*	x:{_,_String},{expressions},*)
(*	x:{{_,_String}..},expressions,*)
(*	_, False];*)
(*If[!expressionpairs,(* treat it as a list of expressions with no strings *)*)
(*expressionpairs=If[Length[names]==1,{expressions,names},*)
(*Transpose[{expressions,names}]]];*)
(*Print[opts/.Expressions->1];*)
(*ExportEOMfuck[{mass,rhs},opts]*)
(*];*)


(* ::Input:: *)
(*ExportEOM[eom_,opts___]:=ExportEOM[MassMatrix[eom],opts];*)
(**)


ExportEOM[{mass_List?MatrixQ,rhs_List?VectorQ},opts___]:=Module[{ostrm,outputfile,quSubs,trigSubs,linearforceterms,forces,forcenames,myopts,zeeusage,expressions,zeeindices,UZ=Identity,symmass=mass,mfilename},
(* write a state-derivative function in Matlab, with mass matrix, right-hand side, and
extras like forces and additional expressions *)
zeeusage=Zees/.{opts}/.Options[ExportEOM];
If[zeeusage===None,UZ=UnZee];
outputfile=OutputFile/.{opts}/.Options[ExportEOM];
mfilename=mFileName/.{opts}/.If[outputfile===None,Options[ExportEOM],mFileName->FileBaseName[outputfile]];(* mFileName if given, else strip OutputFile, else default *)
forces=Forces/.{opts}/.Options[ExportEOM]/.None->{};
(* forces can be a list of forces, or a list of forces and force names, e.g. {{Rx,"Rx"},{Ry,"Ry"}} *)
expressions=Expressions/.{opts}/.Options[ExportEOM]/.None->{};
(* expressions to export are in the form {expression, "name"} where the latter is
the variable name it's given in Matlab, or a list of these pairs *)
expressions=UZ[Replace[expressions,x:{_,_String}->{x}]];(* If called with just {expr, "name"}, wrap it in a list because we always expect a 2-d array *)
If[VectorQ[forces],forcenames=forces, (* else *)
If[MatrixQ[forces],{forces,forcenames}=Transpose[forces]]];
 (* Force mass matrix symmetric *)
Do[symmass[[i,j]]=symmass[[j,i]],{i,Length[symmass]},{j,i}];
zeeindices=zeeusage/.{Minimum->MinimumZee[{symmass,rhs,expressions}],None->{{},{}},All->{All,All}};
ostrm=If[outputfile=!=None, (* there is a file given *)
Print["Opening file \"",outputfile,"\""];
Print["Current directory: ",Directory[]];
Print["Writing function named ", mfilename];
OpenWrite[outputfile], (* else *)
$Output]; (* Okay ready to print out expressions *)
quSubs=formqureplacements;
trigSubs=formtrigreplacements[{symmass,rhs},opts];
myopts=Sequence[qusubs->quSubs,TrigSubs->trigSubs]; (* options generated automatically *)
WriteString[ostrm,ToString[StringForm["function xdot = `1`(t, x)\n\n",mfilename]]]; (* intended for ode45 call *)
(* Print out some nice comment reminders *)
WriteString[ostrm,"% State derivative code generated by Dynamics Workbench " <> DateString[]<>"\n"];
WriteString[ostrm,"% Define constants\n\n"];
(* Print out the forces, if they exist *)
WriteString[ostrm,"% Define forces: "];
(* Convert the force variable into CForm, separate by commas, join these strings,
and put the whole thing in brackets to output *)
If[Length[forces]>0, WriteString[ostrm,"[",
StringJoin[Riffle[Map[ToString[CForm[#]]&,forces],", "]], "]'"]];
WriteString[ostrm,"\n\n% State assignments\n"];
WriteString[ostrm,(* export the state assignments, q1 = x(1), etc. *)
StringJoin[Map[ToString[StringForm["q`1` = x(`1`); ",#]]&,qdofs]],"\n"];WriteString[ostrm,StringJoin[Map[ToString[StringForm["u`1` = x(`2`); ", #,#+Length[qdofs]]]&,udofs]],"\n\n"];
(* export the trig assignments, s1 = sin(q1), etc. *)
WriteString[ostrm,StringReplace[ (* strip out the quote marks *)
StringJoin[Map[ToString[StringForm["`1` = `2`; ", #[[2]],CForm[#[[1]]/.quSubs/.{Sin->"sin",Cos->"cos",Tan->"tan",Sec->"sec",Csc->"csc"}]]]&,trigSubs]],"\""->""],"\n\n"];
(* Begin exporting expressions *)
If[(Length[zees]+Length[zeesC]>0) && zeeusage=!=None&&zeeindices=!={{},{}}, (* Zees on, used, non-empty *)
WriteString[ostrm,"% ZeesC (short-hand constant expressions)\n",ExportExpression[zeesC,"ZC",myopts,opts,Indices->zeeindices[[1]]],"\n\n"];
WriteString[ostrm,"% Zees (short-hand expressions)\n",ExportExpression[zees,"Z",myopts,opts,Indices->zeeindices[[2]]],"\n"];
];
(* export something of the form{MM,"MM"} *)
WriteString[ostrm,"% Mass Matrix\n",ExportExpression[UZ[symmass],"MM",myopts,opts],"\n"];
If[Length[forces]>0, (* if there are listed forces, export the linear terms for them *)
linearforceterms = LinearTerms[rhs,forces,opts];
WriteString[ostrm,"% Force Matrix\n",ExportExpression[UZ[linearforceterms],"FF",myopts,opts],"\n"];
WriteString[ostrm,"% other righthand side terms\n",ExportExpression[UZ[UnLinearTerms[rhs,forces]],"rhs",myopts,opts],"\n"];
WriteString[ostrm,"rhs = rhs + FF*forces;\n"], (* else just give the whole rhs *)
WriteString[ostrm,"% righthand side terms\n",ExportExpression[UZ[rhs],"rhs",myopts,opts],"\n"]];
WriteString[ostrm,"udot = MM\\rhs;\n"];
WriteString[ostrm,ToString[StringForm["xdot = [x(`1`+1:2*`1`); udot];",Length[qdofs]]],"\n"];
(* write out extra expressions, if any *)
If[Length[expressions]>0,WriteString[ostrm,Apply[Sequence,Flatten[Map[{"\n",ExportExpression[#[[1]],#[[2]],myopts,opts],"\n"}&,expressions]]]]];
WriteString[ostrm,"\nend % function\n"];If[outputfile=!=None,Print["Closing file \"",outputfile,"\""];Close[ostrm]];
];
ExportEOM[eom_,opts___]:=ExportEOM[MassMatrix[eom],opts];
ExportCode[{mass_List?MatrixQ,rhs_List?VectorQ},opts___]:=ExportEOM[{mass,rhs},opts];
ExportCode[eom_,opts___]:=ExportEOM[MassMatrix[eom],opts];


(* ::Text:: *)
(*ExportVariables writes a Matlab function to evaluate the specified variables. Features include substitution of q's and u's for the state x, and substitution of shorthand trig functions such as s1 for sin(q(1)). The variables can simply be listed as in ExportVariables[{var1,var2,...}], in which case they will be exported with names as given. Optional text names can also be supplied by pairing with each variable, e.g. ExportVariables[{{var1,"Var1"},{var2,"Var2"}}]. The output of can be sent to a file with OutputFile->"filename", which is automatically produced in the current directory. The Matlab function is named after the OutputFile if specified, with default "exportedvariables". It can also be explicitly set with option mFileName->"name". Formatting options include EndLine->"\n", Continuator->"...", LineLength->78, Deliminator->" ", which specify the end-of-line character, line continuator, maximum line length of the output, and the deliminators where line breaks can occur (default to space character). Zees-> All, None, Minimum*)


(* ::Input:: *)
(*(* GetNames allows the user to just list a bunch of variables without pairing with a text string, and we will try to figure out the name of the Mathematica variable and return it as a string. So inputs like GetNames[{x, y}] will return {{x,"x"},{y,"y"}}.  *)*)
(*SetAttributes[GetNames,{HoldAll,Listable}];*)
(*GetNames[variable_]:={ReleaseHold[variable],ToString[HoldForm[variable]]};*)
(**)
(*SetAttributes[ExportVariables,HoldFirst];*)
(*ExportVariables[variables:{{_,_String}..},opts___]:=ExportDWExpressions[ReleaseHold[variables],opts];*)
(*ExportVariables[variables:{_,_String},opts___]:=ExportDWExpressions[ReleaseHold[variables],opts];*)
(*ExportVariables[variables_,opts___]:=ExportDWExpressions[GetNames[variables],opts];*)


Options[ExportDWExpressions]={OutputFile->None,Zees->Minimum,mFileName->"exportedvariables"};(*ExportDWExpressions[variables:{{_,_String}..},opts___]:=*)
ExportDWExpressions[variables_,opts___]:=Module[{ostrm,outputfile,quSubs,trigSubs,linearforceterms,forces,forcenames,myopts,zeeusage,expressions,zeeindices,UZ=Identity,mfilename,outvarnames},
(* write a function in Matlab, with the variables requested *)
zeeusage=Zees/.{opts}/.Options[ExportDWExpressions];
If[zeeusage===None,UZ=UnZee];
outputfile=OutputFile/.{opts}/.Options[ExportDWExpressions];
mfilename=mFileName/.{opts}/.If[outputfile===None,Options[ExportDWExpressions],mFileName->FileBaseName[outputfile]];(* mFileName if given, else strip OutputFile, else default *)
(* variables\.08 to export are in the form {expression, "name"} where the latter is
the variable name it's given in Matlab, or a list of these pairs *)
expressions=UZ[variables];
zeeindices=zeeusage/.{Minimum->MinimumZee[expressions],None->{{},{}},All->{All,All}};
ostrm=If[outputfile=!=None, (* there is a file given *)
Print["Opening file \"",outputfile,"\""];
Print["Current directory: ",Directory[]];
Print["Writing function named ", mfilename];
OpenWrite[outputfile], (* else *)
$Output]; (* Okay ready to print out expressions *)
quSubs=formqureplacements;
trigSubs=formtrigreplacements[variables,opts];
myopts=Sequence[qusubs->quSubs,TrigSubs->trigSubs]; (* options generated automatically *)
outvarnames = Row[Transpose[variables][[2]],","];
WriteString[ostrm,ToString[StringForm["function [`1`] = `2`(t, x)\n\n",outvarnames,mfilename]]]; (* Print out some nice comment reminders *)
WriteString[ostrm,"% Variables exported by Dynamics Workbench " <> DateString[]<>"\n"];
WriteString[ostrm,"\n\n% State assignments\n"];
WriteString[ostrm,(* export the state assignments, q1 = x(1), etc. *)
StringJoin[Map[ToString[StringForm["q`1` = x(`1`); ",#]]&,qdofs]],"\n"];WriteString[ostrm,StringJoin[Map[ToString[StringForm["u`1` = x(`2`); ", #,#+Length[qdofs]]]&,udofs]],"\n\n"];
(* export the trig assignments, s1 = sin(q1), etc. *)
WriteString[ostrm,StringReplace[ (* strip out the quote marks *)
StringJoin[Map[ToString[StringForm["`1` = `2`; ", #[[2]],CForm[#[[1]]/.quSubs/.{Sin->"sin",Cos->"cos",Tan->"tan",Sec->"sec",Csc->"csc"}]]]&,trigSubs]],"\""->""],"\n\n"];
(* Begin exporting expressions *)
If[(Length[zees]+Length[zeesC]>0) && zeeusage=!=None&&zeeindices=!={{},{}}, (* Zees on, used, non-empty *)
WriteString[ostrm,"% ZeesC (short-hand constant expressions)\n",ExportExpression[zeesC,"ZC",myopts,opts,Indices->zeeindices[[1]]],"\n\n"];
WriteString[ostrm,"% Zees (short-hand expressions)\n",ExportExpression[zees,"Z",myopts,opts,Indices->zeeindices[[2]]],"\n"];
];
(* write out the variables *)
If[Length[expressions]>0,WriteString[ostrm,Apply[Sequence,Flatten[Map[{"\n",ExportExpression[#[[1]],#[[2]],myopts,opts],"\n"}&,expressions]]]]];
WriteString[ostrm,"\nend % function\n"];If[outputfile=!=None,Print["Closing file \"",outputfile,"\""];Close[ostrm]];
];
ExportDWExpressions[variables:{_,_String},opts___]:=ExportDWExpressions[{variables},opts];


(* ::Text:: *)
(*Deprecated ExportCode in here*)


(* ::Input:: *)
(*ExportCode[{mass_List?MatrixQ,rhs_List?VectorQ},opts___]:=Module[{ostrm,outputfile,quSubs,trigSubs,linearforceterms,forces,forcenames,myopts},*)
(*(* write a state-derivative function in Matlab, with mass matrix, right-hand side, and*)
(*extras like forces and additional expressions *)*)
(*outputfile=OutputFile/.{opts}/.Options[ExportEquation];*)
(*forces=Forces/.{opts}/.Options[ExportEquation]/.None->{};*)
(*(* forces can be a list of forces, or a list of forces and force names, e.g. {{Rx,"Rx"},{Ry,"Ry"}} *)*)
(*If[VectorQ[forces],forcenames=forces, (* else *)*)
(*If[MatrixQ[forces],{forces,forcenames}=Transpose[forces]]];*)
(*expressions=Replace[Expressions/.{opts}/.Options[ExportEquation]/.None->{}, x:{_,_String}->{x}];*)
(*ostrm=If[outputfile=!=None, (* there is a file given *)*)
(*Print["Opening file \"",outputfile,"\""];*)
(*Print["Current directory: ",Directory[]];*)
(*OpenWrite[outputfile], (* else *)*)
(*$Output]; (* Okay ready to print out expressions *)*)
(*quSubs=formqureplacements;*)
(*trigSubs=formtrigreplacements[{mass,rhs},opts];*)
(*myopts=Sequence[qusubs->quSubs,TrigSubs->trigSubs]; (* options generated automatically *)*)
(*WriteString[ostrm,"function xdot = fstatederivative(t, x)\n\n"]; (* intended for ode45 call *)*)
(*(* Print out some nice comment reminders *)*)
(*WriteString[ostrm,"% Define constants\n\n"];*)
(*(* Print out the forces, if they exist *)*)
(*WriteString[ostrm,"% Define forces: "];*)
(*(* Convert the force variable into CForm, separate by commas, join these strings,*)
(*and put the whole thing in brackets to output *)*)
(*If[Length[forces]>0, WriteString[ostrm,"[",*)
(*StringJoin[Riffle[Map[ToString[CForm[#]]&, forces],", "]], "]'"]];*)
(*WriteString[ostrm,"\n\n% State assignments\n"];*)
(*WriteString[ostrm,(* export the state assignments, q1 = x(1), etc. *)*)
(*StringJoin[Map[ToString[StringForm["q`1` = x(`1`); ",#]]&,qdofs]],"\n"];WriteString[ostrm,StringJoin[Map[ToString[StringForm["u`1` = x(`2`); ", #,#+Length[qdofs]]]&,udofs]],"\n\n"];*)
(*(* export the trig assignments, s1 = sin(q1), etc. *)*)
(*WriteString[ostrm,StringReplace[ (* strip out the quote marks *)*)
(*StringJoin[Map[ToString[StringForm["`1` = `2`; ", #[[2]],CForm[#[[1]]/.quSubs/.{Sin->"sin",Cos->"cos",Tan->"tan",Sec->"sec",Csc->"csc"}]]]&,trigSubs]],"\""->""],"\n\n"];*)
(*(* pre-allocate mass matrix, force matrix, and otherrhs *)*)
(*If[ZF===Zee, (* Zees on, need to write pre-allocate zees *)*)
(*WriteString[ostrm,ToString[StringForm["ZC = zeros(1,`1`); ",Length[zeesC]]]];WriteString[ostrm,ToString[StringForm["Z = zeros(1,`1`); ",Length[zees]]],"\n"]*)
(*];*)
(*WriteString[ostrm,ToString[StringForm["MM = zeros(`1`,`1`); ",Length[udofs]]]];*)
(*If[Length[forces]>0, (* a force matrix needed *)*)
(*linearforceterms = LinearTerms[rhs,forces,opts];*)
(*WriteString[ostrm,ToString[StringForm["FF = zeros(`1`,`2`); ",Length[udofs],Dimensions[linearforceterms][[2]] ]]]];*)
(*WriteString[ostrm,ToString[StringForm["rhs = zeros(`1`,1);\n\n",Length[udofs]] ]];*)
(*(* Begin exporting expressions *)*)
(*If[ZF===Zee, (* Zees on, need to write pre-allocate zees *)*)
(*WriteString[ostrm,"% ZeesC (short-hand constant expressions)\n",ExportExpression[zeesC,"ZC",myopts,opts],"\n\n"];*)
(*WriteString[ostrm,"% Zees (short-hand expressions)\n",ExportExpression[zees,"Z",myopts,opts],"\n"];*)
(*];*)
(*(* export something of the form{MM,"MM"} *)*)
(*WriteString[ostrm,"% Mass Matrix\n",ExportExpression[mass,"MM",myopts,opts],"\n"];*)
(*If[Length[forces]>0, (* if there are listed forces, export the linear terms for them *)*)
(*WriteString[ostrm,"% Force Matrix\n",ExportExpression[linearforceterms,"FF",myopts,opts],"\n"];*)
(*WriteString[ostrm,"% other righthand side terms\n",ExportExpression[UnLinearTerms[rhs,forces],"rhs",myopts,opts],"\n"];*)
(*WriteString[ostrm,"rhs = rhs + FF*forces;\n"], (* else just give the whole rhs *)*)
(*WriteString[ostrm,"% righthand side terms\n",ExportExpression[rhs,"rhs",myopts,opts],"\n"]];*)
(*WriteString[ostrm,"udot = MM\\rhs;\n"];*)
(*WriteString[ostrm,ToString[StringForm["xdot = [x(`1`+1:2*`1`); udot];",Length[qdofs]]],"\n"];*)
(*(* write out extra expressions, if any *)*)
(*WriteString[ostrm,Apply[Sequence,Flatten[Map[{"\n",ExportExpression[#[[1]],#[[2]],myopts,opts],"\n"}&,expressions]]]];*)
(*If[outputfile=!=None,Print["Closing file \"",outputfile,"\""];Close[ostrm]];*)
(*];*)
(*ExportCode[eom_,opts___]:=ExportCode[MassMatrix[eom],opts];*)


ExportEquations[expr_,opts___]:=Module[{ostrm,outputfile,quSubs,trigSubs,myopts,minzees,expressions,UZ=Identity,zeeusage,zeeindices},
(* write a function in Matlab that calculates the given expression(s), including zees if necessary *)
zeeusage=Zees/.{opts}/.Options[ExportEOM];
If[zeeusage===None,UZ=UnZee];
expressions=UZ[Replace[expr,x:{_,_String}->{x}]];(* If called with just {expr, "name"}, wrap it in a list because we always expect a 2-d array *)
zeeindices=zeeusage/.{Minimum->MinimumZee[expressions],None->{{},{}},All->{All,All}};
outputfile=OutputFile/.{opts}/.Options[ExportEOM];
 (* expressions to export are in the form {expression, "name"} where the latter is the variable name its given in Matlab *)ostrm=If[outputfile=!=None, (* there is a file given *)
Print["Opening file \"",outputfile,"\""];
Print["Current directory: ",Directory[]];
OpenWrite[outputfile], (* else *)
$Output]; (* Okay ready to print out expressions *)
quSubs=formqureplacements;
trigSubs=formtrigreplacements[expressions,opts];
myopts=Sequence[qusubs->quSubs,TrigSubs->trigSubs]; (* options generated automatically *)
WriteString[ostrm,"function output = equations(t, x)\n\n"]; (* a generic function call *)
(* Print out some nice comment reminders *)
WriteString[ostrm,"% Generated by Dynamics Workbench " <> DateString[]<>"\n"];
WriteString[ostrm,"\n\n% State assignments\n"];
WriteString[ostrm,(* export the state assignments, q1 = x(1), etc. *)
StringJoin[Map[ToString[StringForm["q`1` = x(`1`); ",#]]&,qdofs]],"\n"];WriteString[ostrm,StringJoin[Map[ToString[StringForm["u`1` = x(`2`); ", #,#+Length[qdofs]]]&,udofs]],"\n\n"];
(* export the trig assignments, s1 = sin(q1), etc. *)
WriteString[ostrm,StringReplace[ (* strip out the quote marks *)
StringJoin[Map[ToString[StringForm["`1` = `2`; ", #[[2]],CForm[#[[1]]/.quSubs/.{Sin->"sin",Cos->"cos",Tan->"tan",Sec->"sec",Csc->"csc"}]]]&,trigSubs]],"\""->""],"\n\n"];
WriteString[ostrm, "\n"];(* Begin exporting Zees *)
If[(Length[zees]+Length[zeesC]>0) && zeeusage=!=None&&zeeindices=!={{},{}}, (* Zees on, used, not empty *)
WriteString[ostrm,"% ZeesC (short-hand constant expressions)\n",ExportExpression[zeesC,"ZC",myopts,opts,Indices->zeeindices[[1]]],"\n\n"];
WriteString[ostrm,"% Zees (short-hand expressions)\n",ExportExpression[zees,"Z",myopts,opts,Indices->zeeindices[[2]]],"\n"];
];

(* write out the expressions, if any *)
WriteString[ostrm,Apply[Sequence,Flatten[Map[{"\n",ExportExpression[#[[1]],#[[2]],myopts,opts],"\n"}&,expressions]]]];
(* set outputs *)
WriteString[ostrm,"\n% Set outputs\n"];
WriteString[ostrm,StringJoin[ToString[StringForm["output.`1` = `1`;\n", #]]&/@Transpose[expressions][[2]]]];
WriteString[ostrm,"\nend % function\n"];If[outputfile=!=None,Print["Closing file \"",outputfile,"\""];Close[ostrm]];
];


(* ::Section:: *)
(*Building a Model*)


(* ::Subsection:: *)
(*NewModel*)


(* ::Special1:: *)
(*Whenever a new model is constructed, NewModel must be invoked to set up a new frame tree*)


Options[NewModel]={Zees->Off,DotProductList->Off,Simplify->Off};
NewModel[OptionsPattern[]] := Module[{gif},
	(*Begin["DynWkbnch`"];*)
	qdofs={};
	udofs={};
	zees={}; (* short-hand parts of expressions *)
	zeesC = {}; (* Constant values for zees *)
	Clear[ Drawing,Parents, Kids, Tmtx, Kinematics,
	  Inboard, VelCOM, VelJnt, AngVel, AccCOM, AccJnt, AngAcc,
	  Force, Torque, Mass, Inertia, BodyToJnt, InbToJnt, JntToJnt,
	  GInerFrc, GActFrc, Tstar, CML, Constraints, JntAxes, JntType];
	  
	SetAttributes[Drawing,Listable];
	Drawing[ground]={};
	Drawing[All]:=Drawing[Prepend[Bodies,ground]];
	Parents[ground] = {};
	Kids[ground] = {};
	Tmtx[{}] = IdentityMatrix[3];
	Frames = {ground};
	Bodies = {};
	ground[n_?CoordQ] := PV[{1,ground,n}];
	Kinematics = {};
	Nonholonomic = {};
	Mass[(bodies_)?ListQ] := Mass /@ bodies; 
	PosCOM[ ground ] = 0;
	AngVel[ ground ] = PV[{0,ground,1}];
	AngVel[ frme1_, frme2_ ] := AngVel[frme1] - 
		AngVel[frme2] /; ValueQ[AngVel[frme2]];
	VelCOM[ ground ] = 0;
	VelCOM[(bodies_)?ListQ] := 
  Simplify[Plus @@ (Mass /@ bodies*VelCOM /@ bodies)/
    Plus @@ Mass /@ bodies];
	AngAcc[ ground ] = 0;
	AccCOM[ ground ] = 0;
	AccCOM[(bodies_)?ListQ] := Simplify[Plus@@(Mass/@ bodies*AccCOM /@ bodies)/ Plus @@ Mass /@ bodies];
	Inboard[ ground ] = {};
	Dofs[ ground ] = {{},{}};
	ZF = If[OptionValue[Zees]===On,Zee,Identity];
	If[OptionValue[DotProductList]===On, (* want to memorize dot products *)
		DotProduct[frame1_,coord1_,frame2_,coord2_]=. ]; 
	gif = If[OptionValue[Simplify]===On, GInerFrcS, GInerFrcC];
	GInerFrc[x__] := GInerFrc[x] = gif[x]; 
	GActFrc[x__] := GActFrc[x] = GActFrcC[x];
	Tstar[x__] := Tstar[x] = TstarC[x];
    CML[from_,to_]:=Module[{frame1,frame2},{frame1,frame2}=Sort[{from,to}];
      If[OrderedQ[{from,to}],Identity,Transpose][
      CML[frame1,frame2]=CastMtx[frame1,frame2]]]
	(*End[ ];*)
	]


(* ::Subsection:: *)
(*AddFrame*)


(* ::Subsubsection:: *)
(*Options for AddFrame*)


Options[ AddFrame ] = {Axis->{0,0,1}, Qdof->Automatic,
	Axis1->{0,0,1},Axis2->{0,1,0},Axis3->{0,0,1},Qdof1->Automatic,
	Qdof2->Automatic, Qdof3->Automatic, Qdof4->Automatic, 
	Offset->{{1,0,0},{0,1,0},{0,0,1}},  QOffset->0, QOffsets->{0,0}
	};


(* ::Subsubsection:: *)
(*AddFrame for Hinge joint*)


AddFrame[frme_,basefrm_,Hinge, opts___Rule ] := Module[
 {newqdof, z1, z2, z3},
 
	{newqdof} = ReplaceDofs[ {Qdof} 
			/. {opts} /. Options[ AddFrame ], qdofs ];
	qdofs = Union[ qdofs, {newqdof} ];
	Frames = Union[ Append[Frames, frme] ];
	
	q0 = QOffset /. {opts} /. Options[ AddFrame ];
	
	{z1,z2,z3} = ConvertList[Axis /.
				{opts} /. Options[ AddFrame ], basefrm] ;

	(* Now set up the transformation matrix between frames *)
	Tmtx[frme] = {
  {z1^2 (1-Cos[q[newqdof]-q0]) + Cos[q[newqdof]-q0],
   z1 z2 (1-Cos[q[newqdof]-q0]) + z3 Sin[q[newqdof]-q0],
   z1 z3 (1-Cos[q[newqdof]-q0]) - z2 Sin[q[newqdof]-q0]},
  {z1 z2 (1-Cos[q[newqdof]-q0]) - z3 Sin[q[newqdof]-q0],
   z2^2 (1-Cos[q[newqdof]-q0]) + Cos[q[newqdof]-q0],
   z2 z3 (1-Cos[q[newqdof]-q0]) + z1 Sin[q[newqdof]-q0]},
  {z1 z3 (1-Cos[q[newqdof]-q0]) + z2 Sin[q[newqdof]-q0],
   z2 z3 (1-Cos[q[newqdof]-q0]) - z1 Sin[q[newqdof]-q0],
   z3^2 (1-Cos[q[newqdof]-q0]) + Cos[q[newqdof]-q0]}};
   
	Evaluate[frme][n_?CoordQ] := PV[{1,frme,n}];
	Parents[frme] = Append[ Parents[basefrm], frme];
	Kids[frme] = {};
	Kids[basefrm] = Append[ Kids[basefrm], frme];

	{newqdof}
];


(* ::Subsubsection:: *)
(*AddFrame for Fixed/Slider joint*)


AddFrame[frme_,basefrm_,Fixed, opts___Rule] := Module[
	{mtx, z1, z2, z3}, 
	Frames = Union[ Append[Frames, frme] ]; 
	
	(* Now set up the transformation matrix between frames *)
	Tmtx[frme] = ConvertOffset[ Offset, basefrm ] /. {opts} /.
		 Options[AddFrame]; 

	Evaluate[frme][n_?CoordQ] := PV[{1,frme,n}]; 
	Parents[frme] = Append[ Parents[basefrm], frme]; 
	Kids[frme] = {}; 
	Kids[basefrm] = Append[ Kids[basefrm], frme];
	{}  (* No generalized coordinate, so output is null *)
];


AddFrame[frme_,basefrm_,Slider, opts___Rule] := Module[
	{mtx, z1, z2, z3}, 
	Frames = Union[ Append[Frames, frme] ];
	
	(* Now set up the transformation matrix between frames *)
	Tmtx[frme] = ConvertOffset[ Offset, basefrm ] /. {opts} /.
		 Options[AddFrame];

	Evaluate[frme][n_?CoordQ] := PV[{1,frme,n}];
	Parents[frme] = Append[ Parents[basefrm], frme];
	Kids[frme] = {};
	Kids[basefrm] = Append[ Kids[basefrm], frme];
	{}  (* No generalized coordinate, so output is null *)
];


(* ::Subsubsection:: *)
(*AddFrame for U joint*)


AddFrame[frme_,basefrm_,UJoint, opts___Rule ] := Module[
		{z11, z12, z13, z21, z22, z23,
		 newqdofs, q1, q2},

		{z11,z12,z13} = 
				ConvertList[Axis1, basefrm] /. {opts} /. Options[ AddFrame ];
		{z21,z22,z23} =
	  	ConvertList[Axis2, basefrm] /. {opts} /. Options[ AddFrame ];
	  	(* Note above Axis2 should really already be in list form or at least axis # form. If it's in DW vector form, there's no way to specify the axis correctly. *)
	
		newqdofs = {q1, q2} = ReplaceDofs[ {Qdof1, Qdof2}/.
		 	{opts} /. Options[ AddFrame ], qdofs ] ;
		qdofs = Union[ qdofs, newqdofs ];
	
		Frames = Union[ Append[Frames, frme] ];
	
	{q10,q20} = QOffsets /. {opts} /. Options[ AddFrame ];
	
		(* Now set up the transformation matrix between frames *)
	Tmtx[frme] = {
	 {z21^2 (1-Cos[q[q2]-q20]) + Cos[q[q2]-q20],
   z21 z22 (1-Cos[q[q2]-q20]) + z23 Sin[q[q2]-q20],
   z21 z23 (1-Cos[q[q2]-q20]) - z22 Sin[q[q2]-q20]},
  {z21 z22 (1-Cos[q[q2]-q20]) - z23 Sin[q[q2]-q20],
   z22^2 (1-Cos[q[q2]-q20]) + Cos[q[q2]-q20],
   z22 z23 (1-Cos[q[q2]-q20]) + z21 Sin[q[q2]-q20]},
  {z21 z23 (1-Cos[q[q2]-q20]) + z22 Sin[q[q2]-q20],
   z22 z23 (1-Cos[q[q2]-q20]) - z21 Sin[q[q2]-q20],
   z23^2 (1-Cos[q[q2]-q20]) + Cos[q[q2]-q20]}} .{
		{z11^2 (1-Cos[q[q1]-q10]) + Cos[q[q1]-q10],
   z11 z12 (1-Cos[q[q1]-q10]) + z13 Sin[q[q1]-q10],
   z11 z13 (1-Cos[q[q1]-q10]) - z12 Sin[q[q1]-q10]},
  {z11 z12 (1-Cos[q[q1]-q10]) - z13 Sin[q[q1]-q10],
   z12^2 (1-Cos[q[q1]-q10]) + Cos[q[q1]-q10],
   z12 z13 (1-Cos[q[q1]-q10]) + z11 Sin[q[q1]-q10]},
  {z11 z13 (1-Cos[q[q1]-q10]) + z12 Sin[q[q1]-q10],
   z12 z13 (1-Cos[q[q1]-q10]) - z11 Sin[q[q1]-q10],
   z13^2 (1-Cos[q[q1]-q10]) + Cos[q[q1]-q10]}};
	
		Evaluate[frme][n_?CoordQ] := PV[{1,frme,n}];
		Parents[frme] = Append[ Parents[basefrm], frme];
		Kids[frme] = {};
		Kids[basefrm] = Append[ Kids[basefrm], frme];
		newqdofs
];


(* ::Subsubsection:: *)
(*AddFrame for Ball joint*)


AddFrame[frme_,basefrm_,Ball | SixDOF, opts___Rule ] := 
  Module[
	{newqdofs, q1, q2, q3, q4},
	
	newqdofs = {q1,q2,q3,q4} = ReplaceDofs[ 
		{Qdof1,Qdof2,Qdof3,Qdof4} /.
		{opts} /. Options[ AddFrame ], qdofs ];
	qdofs = Union[ qdofs, newqdofs ];
	
	Frames = Union[ Append[Frames, frme] ];

    (* Now set up the transformation matrix between frames *)	
	Tmtx[frme] = {
	  {1 - 2 q[q2]^2 - 2 q[q3]^2,
	   2 (q[q1] q[q2] + q[q3] q[q4]),
	   2 (q[q3] q[q1] - q[q2] q[q4])},
	  {2 (q[q1] q[q2] - q[q3] q[q4]),
	   1 - 2 q[q3]^2 - 2 q[q1]^2,
	   2 (q[q2] q[q3] + q[q1] q[q4])},
	  {2 (q[q3] q[q1] + q[q2] q[q4]),
       2 (q[q2] q[q3] - q[q1] q[q4]),
	   1 - 2 q[q1]^2 - 2 q[q2]^2}};
	
	Evaluate[frme][n_?CoordQ] := PV[{1,frme,n}];
	Parents[frme] = Append[ Parents[basefrm], frme];
	Kids[frme] = {};
	Kids[basefrm] = Append[ Kids[basefrm], frme];
	newqdofs
];


(* ::Subsubsection:: *)
(*AddFrame for Gimbal joint*)


AddFrame[frme_,basefrm_,Gimbal, opts___Rule ] := Module[
		{z11, z12, z13, z21, z22, z23, z31, z32, z33,
		 newqdofs, q1, q2, q3},

		{z11,z12,z13} = 
				ConvertList[Axis1, basefrm] /. {opts} /. Options[ AddFrame ];
		{z21,z22,z23} =
	  	ConvertList[Axis2, basefrm] /. {opts} /. Options[ AddFrame ];
		{z31,z32,z33} =
				ConvertList[Axis3, basefrm] /. {opts} /. Options[ AddFrame ];
	
		newqdofs = {q1, q2, q3} = ReplaceDofs[ {Qdof1, Qdof2,
		  Qdof3 } /.	{opts} /. Options[ AddFrame ], qdofs ] ;
		qdofs = Union[ qdofs, newqdofs ];
	
		Frames = Union[ Append[Frames, frme] ];
	
		(* Now set up the transformation matrix between frames *)
	Tmtx[frme] = {
	 {z31^2 (1-Cos[q[q3]]) + Cos[q[q3]],
   z31 z32 (1-Cos[q[q3]]) + z33 Sin[q[q3]],
   z31 z33 (1-Cos[q[q3]]) - z32 Sin[q[q3]]},
  {z31 z32 (1-Cos[q[q3]]) - z33 Sin[q[q3]],
   z32^2 (1-Cos[q[q3]]) + Cos[q[q3]],
   z32 z33 (1-Cos[q[q3]]) + z31 Sin[q[q3]]},
  {z31 z33 (1-Cos[q[q3]]) + z32 Sin[q[q3]],
   z32 z33 (1-Cos[q[q3]]) - z31 Sin[q[q3]],
   z33^2 (1-Cos[q[q3]]) + Cos[q[q3]]}} . {
  {z21^2 (1-Cos[q[q2]]) + Cos[q[q2]],
   z21 z22 (1-Cos[q[q2]]) + z23 Sin[q[q2]],
   z21 z23 (1-Cos[q[q2]]) - z22 Sin[q[q2]]},
  {z21 z22 (1-Cos[q[q2]]) - z23 Sin[q[q2]],
   z22^2 (1-Cos[q[q2]]) + Cos[q[q2]],
   z22 z23 (1-Cos[q[q2]]) + z21 Sin[q[q2]]},
  {z21 z23 (1-Cos[q[q2]]) + z22 Sin[q[q2]],
   z22 z23 (1-Cos[q[q2]]) - z21 Sin[q[q2]],
   z23^2 (1-Cos[q[q2]]) + Cos[q[q2]]}} . {
  {z11^2 (1-Cos[q[q1]]) + Cos[q[q1]],
   z11 z12 (1-Cos[q[q1]]) + z13 Sin[q[q1]],
   z11 z13 (1-Cos[q[q1]]) - z12 Sin[q[q1]]},
  {z11 z12 (1-Cos[q[q1]]) - z13 Sin[q[q1]],
   z12^2 (1-Cos[q[q1]]) + Cos[q[q1]],
   z12 z13 (1-Cos[q[q1]]) + z11 Sin[q[q1]]},
  {z11 z13 (1-Cos[q[q1]]) + z12 Sin[q[q1]],
   z12 z13 (1-Cos[q[q1]]) - z11 Sin[q[q1]],
   z13^2 (1-Cos[q[q1]]) + Cos[q[q1]]}};
	
		Evaluate[frme][n_?CoordQ] := PV[{1,frme,n}];
		Parents[frme] = Append[ Parents[basefrm], frme];
		Kids[frme] = {};
		Kids[basefrm] = Append[ Kids[basefrm], frme];
		newqdofs
];


(* ::Subsection:: *)
(*AddBody*)


(* ::Subsubsection:: *)
(*Options for AddBody*)


Options[ AddBody ] = {Axis->{0,0,1}, Qdof->Automatic,
  Udof->Automatic, Udof1->Automatic, Udof2->Automatic,
  Udof3->Automatic,Udof4->Automatic,Udof5->Automatic,
  Udof6->Automatic,
  Qdof1->Automatic, Qdof2->Automatic,Qdof3->Automatic,
  Qdof4->Automatic, Qdof5->Automatic, Qdof6->Automatic,
  Qdof7->Automatic,
  Axis1->{0,0,1}, Axis2->{0,1,0},Axis3->{0,0,1},
  Taxis->{1,0,0},
  TAxis1->{1,0,0},TAxis2->{0,1,0},TAxis3->{0,0,1},
	BodyToJnt->0, InbToJnt->0, JntToJnt->0,Mass->0, Inertia->0,
	Frme->Automatic, Basefrm->Automatic,
	RelativeTo->Automatic, Drawing->Automatic};


(* ::Subsubsection:: *)
(*AddBody for Hinge joint*)


AddBody[ body_, inboard_, Hinge, opts___Rule ] := Module[ 
		{BtJ, ItJ, newqdof, newudof, frme,
  	basefrm, axis, com, jnt, drwing},

	 frme = Frme /. {opts} /. Options[ AddBody ];
		If[ frme == Automatic, frme = body ];	
		basefrm = RelativeTo /. {opts} /. Options[ AddBody ];
		If[ basefrm == Automatic, basefrm = inboard ];
  {newqdof} = AddFrame[ frme, basefrm, Hinge, opts ];
	 {newudof} = ReplaceDofs[{Udof} /. {opts} /.
	   Options[ AddBody ], udofs ];
	 udofs = Union[ udofs, {newudof} ];

		axis = ConvertVec[Axis, basefrm] /. 
				{opts} /. Options[ AddFrame ];
	JntAxes[body] = axis;
	JntType[body] = Hinge;
	
 	(* Handle the bodies *)
	
 	AppendTo[ Bodies, body ];
 	Inboard[body] = inboard;
 	{BtJ, ItJ, Mass[body], Inertia[body]} =
				{BodyToJnt, InbToJnt, Mass, ConvertInertia[Inertia,frme]} /.
				{opts} /. Options[ AddBody ];
		BodyToJnt[body] = BtJ;
		InbToJnt[body] = ItJ;
		JntToJnt[body] = 0;
	
 	AppendTo[ Kinematics, u[newudof] == qt[newqdof]'[t] ];
 	Kinematics = Union[Kinematics];
 	AngVel[ body ] = AngVel[ basefrm ] + 
 	  u[newudof] * axis;
 	VelJnt[ body ] = VelCOM[ inboard ] +
 	  Cross[ AngVel[inboard], ItJ ];
 	VelCOM[ body ] = VelJnt[ body ] +
 	  Cross[ AngVel[body], -BtJ ] //ZF;
 	AngAcc[ body ] = Dt[ AngVel[body], t];
 	AccJnt[ body ] = AccCOM[ inboard ] +
 	  Cross[ AngVel[inboard], Cross[ AngVel[inboard], ItJ ] ] +
 	  Cross[ AngAcc[inboard], ItJ ];
 	AccCOM[ body ] = AccJnt[ body ] + 
 	  Cross[ AngVel[body], Cross[ AngVel[body], -BtJ ] ] +
 	  Cross[ AngAcc[body], -BtJ ] //ZF;
 	Force[ body ] = {};
 	Torque[ body ] = {};
 	com = Thread[PosCOM[ body ].{ground[1],ground[2],ground[3]}]; (* Cartesian coords of com *)
 	jnt = com+Thread[BtJ.{ground[1],ground[2],ground[3]}]; (* Cartesian coords of joint *)
 	drwing = Drawing /.{opts}/.Options[ AddBody ]; (* if a Drawing option specified, use it *)
 	Drawing[ body ] = If[ drwing === Automatic,
 	  If[ inboard != ground, (* then add a line from inboard COM to the joint *)
 	    AppendTo[Drawing[inboard],Line[Thread[PosCOM[inboard].{ground[1],ground[2],ground[3]}],jnt]]
 	  ];
 	  {Point[jnt],Line[{jnt,com}],Point[com]}
 	  , (* else use the optional drawing, centered around the body's com *)
 	  Translate[drwing,com]
 	];
 	Dofs[ body ] = {{newqdof}, {newudof}}
];


(* ::Subsubsection:: *)
(*AddBody for Fixed joint*)


AddBody[ body_, inboard_, Fixed, opts___Rule ] := Module[ 
	 {BtJ, ItJ, frme,
  	basefrm, taxis, com, jnt, drwing},
		
		frme = Frme /. {opts} /. Options[ AddBody ];
		If[ frme == Automatic, frme = body ];	
		basefrm = RelativeTo /. {opts} /. Options[ AddBody ];
		If[ basefrm == Automatic, basefrm = inboard ];
  AddFrame[ frme, basefrm, Fixed, opts ];
  JntAxes[body] = {};
  JntType[body] = Fixed;
	
 	(* Handle the bodies *)
	
 	AppendTo[ Bodies, body ];
 	Inboard[body] = inboard;
  BtJ = BodyToJnt /. {opts} /. Options[ AddBody ];
 	ItJ = InbToJnt /. {opts} /. Options[ AddBody ];
		BodyToJnt[body] = BtJ;
		InbToJnt[body] = ItJ;
		JntToJnt[body] = 0;
 	Mass[body] = Mass /. {opts} /. Options[ AddBody ];
 	Inertia[body] = ConvertInertia[Inertia,frme] /. {opts} /. Options[ AddBody ];
	
 	AngVel[ body ] = AngVel[ inboard ];
 	VelJnt[ body ] = VelCOM[ inboard ] +
 	  Cross[ AngVel[inboard], ItJ ];
 	VelCOM[ body ] = VelJnt[ body ] +
 	  Cross[ AngVel[body], -BtJ ] //ZF;
 	AngAcc[ body ] = Dt[ AngVel[body], t];
 	AccJnt[ body ] = AccCOM[ inboard ] +
 	  Cross[ AngVel[inboard], Cross[ AngVel[inboard], ItJ ] ] +
 	  Cross[ AngAcc[inboard], ItJ ];
 	AccCOM[ body ] = AccJnt[ body ] + 
 	  Cross[ AngVel[body], Cross[ AngVel[body], -BtJ ] ] +
 	  Cross[ AngAcc[body], -BtJ ] //ZF;
 	Force[ body ] = {};
 	Torque[ body ] = {};
 	com = Thread[PosCOM[ body ].{ground[1],ground[2],ground[3]}]; (* Cartesian coords of com *)
 	jnt = com+Thread[BtJ.{ground[1],ground[2],ground[3]}]; (* Cartesian coords of joint *)
    drwing = Drawing /.{opts}/.Options[ AddBody ]; (* if a Drawing option specified, use it *)
 	Drawing[ body ] = If[ drwing === Automatic,
 	  If[ inboard != ground, (* then add a line from inboard COM to the joint *)
 	    AppendTo[Drawing[inboard],Line[Thread[PosCOM[inboard].{ground[1],ground[2],ground[3]}],jnt]]
 	  ];
 	  {Point[jnt],Line[{jnt,com}],Point[com]}
 	  , (* else use the optional drawing, centered around the body's com *)
 	  Translate[drwing,com]
 	];
Dofs[ body ] = {{}, {}}
];


(* ::Subsubsection:: *)
(*AddBody for Slider joint*)


AddBody[ body_, inboard_, Slider, opts___Rule ] := Module[ 
	 {BtJ, ItJ, newqdof, newudof, frme,
  	basefrm, taxis, drwing,com,jnt},
		
		frme = Frme /. {opts} /. Options[ AddBody ];
		If[ frme == Automatic, frme = body ];	
		basefrm = RelativeTo /. {opts} /. Options[ AddBody ];
		If[ basefrm == Automatic, basefrm = inboard ];
  {newqdof} = ReplaceDofs[{Qdof} /. {opts} /.
  	 Options[ AddBody ], qdofs];
  qdofs = Union[ qdofs, {newqdof} ]; 
  AddFrame[ frme, basefrm, Fixed, opts ]; 
	 {newudof} = ReplaceDofs[{Udof} /. 
	 		{opts} /. Options[ AddBody ], udofs];
	 udofs = Union[ udofs, {newudof} ]; 
	
		taxis = ConvertVec[TAxis, basefrm] /. {opts} /.
			 Options[ AddFrame ]; 
	JntAxes[body] = taxis;
	JntType[body] = Slider;

 	(* Handle the bodies *)
	
 	AppendTo[ Bodies, body ]; 
 	Inboard[body] = inboard; 
  BtJ = BodyToJnt /. {opts} /. Options[ AddBody ];
 	ItJ = InbToJnt /. {opts} /. Options[ AddBody ];
		BodyToJnt[body] = BtJ;
		InbToJnt[body] = ItJ;
		JntToJnt[body] = q[newqdof] taxis;
 	Mass[body] = Mass /. {opts} /. Options[ AddBody ];
 	Inertia[body] = ConvertInertia[Inertia,frme] /. {opts} /. Options[ AddBody ];
	
 	AppendTo[ Kinematics, u[newudof] == qt[newqdof]'[t] ];
 	Kinematics = Union[Kinematics];
 	AngVel[ body ] = AngVel[ inboard ]; 
 	VelJnt[ body ] = VelCOM[ inboard ] +
 	  Cross[ AngVel[inboard], ItJ ];
 	VelJnt2[ body ] = VelJnt[ body ] + Cross[ AngVel[ body ],
 	  q[newqdof] taxis ] + u[newudof] taxis;
 	(*VelCOM[ body ] = VelJnt2[ body ] +
 	  Cross[ AngVel[body], -BtJ ]//ZF;*)
 	VelCOM[ body ] = VelJnt[ body ] + Cross[ AngVel[body], q[newqdof] taxis - BtJ ] + u[newudof] taxis //ZF;
 	AngAcc[ body ] = Dt[ AngVel[body], t];
 	AccJnt[ body ] = AccCOM[ inboard ] +
 	  Cross[ AngVel[inboard], Cross[ AngVel[inboard], ItJ ] ] +
 	  Cross[ AngAcc[inboard], ItJ ];
 	(*AccJnt2[ body ] = AccJnt[ body] +
 		 Cross[ AngVel[inboard], Cross[ AngVel[inboard], 
 		   q[newqdof] taxis ] ] +
 		 Cross[ AngAcc[inboard], q[newqdof] taxis ] + u[newudof]' taxis + 
 		 2 Cross[ AngVel[body], u[newudof] taxis ];
 	AccCOM[ body ] = AccJnt2[ body ] + 
 	  Cross[ AngVel[body], Cross[ AngVel[body], -BtJ ] ] +
 	  Cross[ AngAcc[body], -BtJ ]//ZF;*)
 	AccCOM[ body ] = AccJnt[ body ] + 
 	  Cross[ AngVel[body], Cross[ AngVel[body], q[newqdof] taxis - BtJ ] ] + 
 	  Cross[ AngAcc[body], q[newqdof] taxis - BtJ ] + u[newudof]' taxis + 
 	  2 Cross[ AngVel[body], u[newudof] taxis ] //ZF;
 	Force[ body ] = {};
 	Torque[ body ] = {}; 
 	com = Thread[PosCOM[ body ].{ground[1],ground[2],ground[3]}]; (* Cartesian coords of com *)
 	jnt = com+Thread[BtJ.{ground[1],ground[2],ground[3]}]; (* Cartesian coords of joint *)
    drwing = Drawing /.{opts}/.Options[ AddBody ]; (* if a Drawing option specified, use it *)
 	Drawing[ body ] = If[ drwing === Automatic,
 	  If[ inboard != ground, (* then add a line from inboard COM to the joint *)
 	    AppendTo[Drawing[inboard],Line[Thread[PosCOM[inboard].{ground[1],ground[2],ground[3]}],jnt]]
 	  ];
 	  {Point[jnt],Line[{jnt,com}],Point[com]}
 	  , (* else use the optional drawing, centered around the body's com *)
 	  Translate[drwing,com]
 	];Dofs[ body ] = {{newqdof}, {newudof}}
];


(* ::Subsubsection:: *)
(*AddBody for UJoint*)


AddBody[ body_, inboard_, UJoint, opts___Rule ] := Module[ 
	 {BtJ, ItJ, newqdofs, newudofs, q1, q2, u1, u2, frme,
   basefrm, axis1, axis2, com, jnt, drwing},
  	     	 
 	frme = Frme /. {opts} /. Options[ AddBody ];
 	If[ frme == Automatic, frme = body ];	
 	basefrm = RelativeTo /. {opts} /. Options[ AddBody ];
 	If[ basefrm == Automatic, basefrm = inboard ];
  newqdofs = {q1, q2} = AddFrame[ frme, basefrm, UJoint, opts ];
  newudofs = {u1, u2} = ReplaceDofs[{Udof1, Udof2} /. 
  		{opts} /. Options[ AddBody ], udofs];
  udofs = Union[ udofs, newudofs ];
	
		axis1 = ConvertVec[Axis1, basefrm] /. {opts} /.
				Options[ AddFrame ];
		axis2 = ConvertVec[Axis2, frme] /. {opts} /.
		  Options[ AddFrame ];
JntAxes[body] = {axis1, axis2};	  
	JntType[body] = UJoint;
	  
 	(* Handle the bodies *)
	
 	AppendTo[ Bodies, body ];
 	Inboard[body] = inboard;
  BtJ = BodyToJnt /. {opts} /. Options[ AddBody ];
 	ItJ = InbToJnt /. {opts} /. Options[ AddBody ];
		BodyToJnt[body] = BtJ;
		InbToJnt[body] = ItJ;
		JntToJnt[body] = 0;
 	Mass[body] = Mass /. {opts} /. Options[ AddBody ];
 	Inertia[body] = ConvertInertia[Inertia,frme] /. {opts} /. Options[ AddBody ];
	
 	AppendTo[ Kinematics, u[u1] == qt[q1]'[t] ];
 	AppendTo[ Kinematics, u[u2] == qt[q2]'[t] ];
 	Kinematics = Union[Kinematics];
 	AngVel[ body ] = AngVel[ inboard ] + 
 	  u[u1] axis1 +
 	  u[u2] axis2;
 	VelJnt[ body ] = VelCOM[ inboard ] +
 	  Cross[ AngVel[inboard], ItJ ];
 	VelCOM[ body ] = VelJnt[ body ] +
 	  Cross[ AngVel[body], -BtJ ]//ZF;
 	AngAcc[ body ] = Dt[ AngVel[body], t];
 	AccJnt[ body ] = AccCOM[ inboard ] +
 	  Cross[ AngVel[inboard], Cross[ AngVel[inboard], ItJ ] ] +
 	  Cross[ AngAcc[inboard], ItJ ];
 	AccCOM[ body ] = AccJnt[ body ] + 
 	  Cross[ AngVel[body], Cross[ AngVel[body], -BtJ ] ] +
 	  Cross[ AngAcc[body], -BtJ ]//ZF;
 	Force[ body ] = {};
 	Torque[ body ] = {};
 	com = Thread[PosCOM[ body ].{ground[1],ground[2],ground[3]}]; (* Cartesian coords of com *)
 	jnt = com+Thread[BtJ.{ground[1],ground[2],ground[3]}]; (* Cartesian coords of joint *)
    drwing = Drawing /.{opts}/.Options[ AddBody ]; (* if a Drawing option specified, use it *)
 	Drawing[ body ] = If[ drwing === Automatic,
 	  If[ inboard != ground, (* then add a line from inboard COM to the joint *)
 	    AppendTo[Drawing[inboard],Line[Thread[PosCOM[inboard].{ground[1],ground[2],ground[3]}],jnt]]
 	  ];
 	  {Point[jnt],Line[{jnt,com}],Point[com]}
 	  , (* else use the optional drawing, centered around the body's com *)
 	  Translate[drwing,com]
 	];
 	Dofs[ body ] = {newqdofs, newudofs}
];


(* ::Subsubsection:: *)
(*AddBody for Gimbal joint*)


AddBody[ body_, inboard_, Gimbal, opts___Rule ] := Module[ 
	 {BtJ, ItJ, newqdofs, newudofs, q1, q2, q3, u1, u2, u3,
	  frme, basefrm, axis1, axis2, axis3, com, jnt, drwing},
  	     	 
 	frme = Frme /. {opts} /. Options[ AddBody ];
 	If[ frme == Automatic, frme = body ];	
 	basefrm = RelativeTo /. {opts} /. Options[ AddBody ];
 	If[ basefrm == Automatic, basefrm = inboard ];
  newqdofs = {q1, q2, q3} = 
  	AddFrame[ frme, basefrm, Gimbal, opts ];
  newudofs = {u1, u2, u3} = 
  	ReplaceDofs[{Udof1, Udof2, Udof3} /. 
  	{opts} /. Options[ AddBody ], udofs];
  udofs = Union[ udofs, newudofs ];
	
	axis1 = ConvertVec[Axis1, basefrm] /. {opts} /.
		Options[ AddFrame ];
	axis2 = ConvertVec[Axis2, basefrm] /. {opts} /.
		Options[ AddFrame ];
	axis3 = ConvertVec[Axis3, basefrm] /. {opts} /.
		Options[ AddFrame ];
	JntAxes[body] = {axis1, axis2, axis3};
	JntType[body] = Gimbal;
	  
 	(* Handle the bodies *)
	
 	AppendTo[ Bodies, body ];
 	Inboard[body] = inboard;
  BtJ = BodyToJnt /. {opts} /. Options[ AddBody ];
 	ItJ = InbToJnt /. {opts} /. Options[ AddBody ];
		BodyToJnt[body] = BtJ;
		InbToJnt[body] = ItJ;
		JntToJnt[body] = 0;
 	Mass[body] = Mass /. {opts} /. Options[ AddBody ];
 	Inertia[body] = ConvertInertia[Inertia,frme] /. {opts} /. Options[ AddBody ];
	
 	AppendTo[ Kinematics, u[u1] == qt[q1]'[t] ];
 	AppendTo[ Kinematics, u[u2] == qt[q2]'[t] ];
 	AppendTo[ Kinematics, u[u3] == qt[q3]'[t] ];
 	Kinematics = Union[Kinematics];
 	AngVel[ body ] = AngVel[ inboard ] + 
 	  u[u1] axis1 +
 	  u[u2] axis2 + u[u3] axis3;
 	VelJnt[ body ] = VelCOM[ inboard ] +
 	  Cross[ AngVel[inboard], ItJ ];
 	VelCOM[ body ] = VelJnt[ body ] +
 	  Cross[ AngVel[body], -BtJ ]//ZF;
 	AngAcc[ body ] = Dt[ AngVel[body], t];
 	AccJnt[ body ] = AccCOM[ inboard ] +
 	  Cross[ AngVel[inboard], Cross[ AngVel[inboard], ItJ ] ] +
 	  Cross[ AngAcc[inboard], ItJ ];
 	AccCOM[ body ] = AccJnt[ body ] + 
 	  Cross[ AngVel[body], Cross[ AngVel[body], -BtJ ] ] +
 	  Cross[ AngAcc[body], -BtJ ]//ZF;
 	Force[ body ] = {};
 	Torque[ body ] = {};
 	 	com = Thread[PosCOM[ body ].{ground[1],ground[2],ground[3]}]; (* Cartesian coords of com *)
 	jnt = com+Thread[BtJ.{ground[1],ground[2],ground[3]}]; (* Cartesian coords of joint *)
    drwing = Drawing /.{opts}/.Options[ AddBody ]; (* if a Drawing option specified, use it *)
 	Drawing[ body ] = If[ drwing === Automatic,
 	  If[ inboard != ground, (* then add a line from inboard COM to the joint *)
 	    AppendTo[Drawing[inboard],Line[Thread[PosCOM[inboard].{ground[1],ground[2],ground[3]}],jnt]]
 	  ];
 	  {Point[jnt],Line[{jnt,com}],Point[com]}
 	  , (* else use the optional drawing, centered around the body's com *)
 	  Translate[drwing,com]
 	];
 	Dofs[ body ] = {newqdofs, newudofs}
];


(* ::Subsubsection:: *)
(*AddBody for Ball Joint*)


AddBody[ body_, inboard_, Ball, opts___Rule ] := Module[ 
	{BtJ, ItJ, q1, q2, q3, q4, u1, u2, u3, frme,
  	basefrm, newqdofs, newudofs,Euler, com, jnt, drwing},
  	 
	frme = Frme /. {opts} /. Options[ AddBody ];
	If[ frme == Automatic, frme = body ];	
	basefrm = RelativeTo /. {opts} /. Options[ AddBody ];
	If[ basefrm == Automatic, basefrm = inboard ];
  newqdofs = {q1,q2,q3,q4} = 
  	AddFrame[ frme, basefrm, Ball, opts ];
  newudofs = {u1, u2, u3} = ReplaceDofs[
  	{Udof1, Udof2, Udof3} /.
  	{opts} /. Options[ AddBody ], udofs];
	udofs = Union[ udofs, newudofs ];
  JntAxes[body] = IdentityMatrix[3];
  JntType[body] = Ball;
	
 	(* Handle the bodies *)
	
 	AppendTo[ Bodies, body ];
 	Inboard[body] = inboard;
  BtJ = BodyToJnt /. {opts} /. Options[ AddBody ];
 	ItJ = InbToJnt /. {opts} /. Options[ AddBody ];
		BodyToJnt[body] = BtJ;
		InbToJnt[body] = ItJ;
		JntToJnt[body] = 0;
 	Mass[body] = Mass /. {opts} /. Options[ AddBody ];
 	Inertia[body] = ConvertInertia[Inertia,frme] /. {opts} /. Options[ AddBody ];
	
	Euler = {{ q[q4], -q[q3],  q[q2]},
		   { q[q3],  q[q4], -q[q1]},
		   {-q[q2],  q[q1],  q[q4]},
		   {-q[q1], -q[q2], -q[q3]}};
		   
 	AppendTo[ Kinematics,Thread[ 
 		{qt[q1]'[t], qt[q2]'[t], qt[q3]'[t], qt[q4]'[t]} ==
 		1/2 Euler . {u[u1], u[u2], u[u3]},List] ];
 	Kinematics = Union[Flatten[Kinematics]];
 	AngVel[ body ] = AngVel[ inboard ] + 
 		PV[{u[u1],frme,1},{u[u2],frme,2},
 		{u[u3],frme,3}];
 	VelJnt[ body ] = VelCOM[ inboard ] +
 	  Cross[ AngVel[inboard], ItJ ];
 	VelCOM[ body ] = VelJnt[ body ] +
 	  Cross[ AngVel[body], -BtJ ]//ZF;
 	AngAcc[ body ] = Dt[ AngVel[body], t];
 	AccJnt[ body ] = AccCOM[ inboard ] +
 	  Cross[ AngVel[inboard], Cross[ AngVel[inboard], ItJ ] ] +
 	  Cross[ AngAcc[inboard], ItJ ];
 	AccCOM[ body ] = AccJnt[ body ] + 
 	  Cross[ AngVel[body], Cross[ AngVel[body], -BtJ ] ] +
 	  Cross[ AngAcc[body], -BtJ ]//ZF;
 	Force[ body ] = {};
 	Torque[ body ] = {};
 	 	com = Thread[PosCOM[ body ].{ground[1],ground[2],ground[3]}]; (* Cartesian coords of com *)
 	jnt = com+Thread[BtJ.{ground[1],ground[2],ground[3]}]; (* Cartesian coords of joint *)
    drwing = Drawing /.{opts}/.Options[ AddBody ]; (* if a Drawing option specified, use it *)
 	Drawing[ body ] = If[ drwing === Automatic,
 	  If[ inboard != ground, (* then add a line from inboard COM to the joint *)
 	    AppendTo[Drawing[inboard],Line[Thread[PosCOM[inboard].{ground[1],ground[2],ground[3]}],jnt]]
 	  ];
 	  {Point[jnt],Line[{jnt,com}],Point[com]}
 	  , (* else use the optional drawing, centered around the body's com *)
 	  Translate[drwing,com]
 	];
 	Dofs[ body ] = {newqdofs, newudofs}
];



(* ::Subsubsection:: *)
(*AddBody for SixDOF Joint*)


AddBody[ body_, inboard_, SixDOF, opts___Rule ] := Module[ 
	{BtJ, ItJ, q1, q2, q3, q4, q5, q6, q7,
	 u1, u2, u3, u4, u5, u6, frme,
  	basefrm, newqdofs, newudofs, taxis1, taxis2, taxis3, Euler, com, jnt, drwing},

	newqdofs = {q1,q2,q3,q4,q5,q6,q7} = ReplaceDofs[
		{Qdof1,Qdof2,Qdof3,Qdof4,Qdof5,Qdof6,Qdof7}
		/. {opts} /. Options[ AddBody ], qdofs];
		
	If[ Length[#] > 0, Message[AddBody::repeatDOF, #] ] & [
  	Intersection[newqdofs, qdofs] ];

	qdofs = Union[ qdofs, newqdofs ]; 	
    	 
	frme = Frme /. {opts} /. Options[ AddBody ];
	If[ frme == Automatic, frme = body ];	
	basefrm = RelativeTo /. {opts} /. Options[ AddBody ];
	If[ basefrm == Automatic, basefrm = inboard ];
  AddFrame[ frme, basefrm, Ball, Qdof1->q4,
  	Qdof2->q5, Qdof3->q6, Qdof4->q7, opts ];
  	
  {taxis1, taxis2, taxis3} = ConvertVec[#,basefrm]& /@
  	{TAxis1, TAxis2, TAxis3} /.
  	{opts} /. Options[ AddBody ];
  	
  newudofs = {u1, u2, u3, u4, u5, u6} = ReplaceDofs[
    {Udof1, Udof2, Udof3, Udof4, Udof5, Udof6} /.
  	{opts} /. Options[ AddBody ], udofs];
  udofs = Union[ udofs, newudofs ];
  JntAxes[body] = IdentityMatrix[3];
	JntType[body] = SixDOF;

 	(* Handle the bodies *)
	
 	AppendTo[ Bodies, body ];
 	Inboard[body] = inboard;
  BtJ = BodyToJnt /. {opts} /. Options[ AddBody ];
 	ItJ = InbToJnt /. {opts} /. Options[ AddBody ];
	BodyToJnt[body] = BtJ;
	InbToJnt[body] = ItJ;
	JntToJnt[body] = 0;
 	Mass[body] = Mass /. {opts} /. Options[ AddBody ];
 	Inertia[body] = ConvertInertia[Inertia,frme] /. {opts} /. Options[ AddBody ];
	
	Euler = {{ q[q7], -q[q6],  q[q5]},
		   { q[q6],  q[q7], -q[q4]},
		   {-q[q5],  q[q4],  q[q7]},
		   {-q[q4], -q[q5], -q[q6]}};
		   
 	AppendTo[ Kinematics, {u[u1],u[u2],u[u3]}=={qt[q1]'[t],qt[q2]'[t],qt[q3]'[t]}];
 	AppendTo[ Kinematics, 
 		{qt[q4]'[t], qt[q5]'[t], qt[q6]'[t], qt[q7]'[t]} ==
 		1/2 Euler . {u[u4], u[u5], u[u6]} ];
 	Kinematics = Union[Kinematics];
 	AngVel[ body ] = AngVel[ inboard ] + 
 		PV[{u[u4],frme,1},{u[u5],frme,2},
 		{u[u6],frme,3}];
 	VelJnt[ body ] = VelCOM[ inboard ] +
 	  Cross[ AngVel[inboard], ItJ ];
 	VelJnt2[ body ] = VelJnt[ body ] + Cross[ AngVel[ inboard ],
 	  q[q1] taxis1 + q[q2] taxis2 + q[q3] taxis3 ] +
 	  u[u1] taxis1 + u[u2] taxis2 + u[u3] taxis3;
 	VelCOM[ body ] = VelJnt2[ body ] +
 	  Cross[ AngVel[body], -BtJ ]//ZF;
 	AngAcc[ body ] = Dt[ AngVel[body], t];
 	AccJnt[ body ] = AccCOM[ inboard ] +
 	  Cross[ AngVel[inboard], Cross[ AngVel[inboard], ItJ ] ] +
 	  Cross[ AngAcc[inboard], ItJ ];
 	AccJnt2[ body ] = AccJnt[ body] +
 		Cross[ AngVel[inboard], Cross[ AngVel[inboard], 
      q[q1] taxis1 + q[q2] taxis2 + q[q3] taxis3 ] ] +
 		Cross[ AngAcc[inboard], 
 			q[q1] taxis1 + q[q2] taxis2 + q[q3] taxis3 ] +
 	 	u[u1]' taxis1 + u[u2]' taxis2 + u[u3]' taxis3 +
 		2 Cross[ AngVel[inboard], 
 		   u[u1] taxis1 + u[u2] taxis2 + u[u3] taxis3 ];
 	AccCOM[ body ] = AccJnt2[ body ] + 
 	  Cross[ AngVel[body], Cross[ AngVel[body], -BtJ ] ] +
 	  Cross[ AngAcc[body], -BtJ ]//ZF;
 	Force[ body ] = {};
 	Torque[ body ] = {};
 	com = Thread[PosCOM[ body ].{ground[1],ground[2],ground[3]}]; (* Cartesian coords of com *)
 	jnt = com+Thread[BtJ.{ground[1],ground[2],ground[3]}]; (* Cartesian coords of joint *)
    drwing = Drawing /.{opts}/.Options[ AddBody ]; (* if a Drawing option specified, use it *)
 	Drawing[ body ] = If[ drwing === Automatic,
 	  If[ inboard != ground, (* then add a line from inboard COM to the joint *)
 	    AppendTo[Drawing[inboard],Line[Thread[PosCOM[inboard].{ground[1],ground[2],ground[3]}],jnt]]
 	  ];
 	  {Point[jnt],Line[{jnt,com}],Point[com]}
 	  , (* else use the optional drawing, centered around the body's com *)
 	  Translate[drwing,com]
 	];
 	Dofs[ body ] = {newqdofs, newudofs}
];



CloseLoop[ body_, inboard_, Hinge, opts___Rule ] := Module[ 
		{BtJ, ItJ, axis, velloopjnt1, velloopjnt2,
		joint, vector1, vector2, angvel, constraints},

	axis = ConvertVec[Axis, body] /. 
		{opts} /. Options[ AddFrame ];
	
 	{BtJ, ItJ} = {BodyToJnt, InbToJnt} /.
		{opts} /. Options[ AddBody ];

	Loop[body] = {inboard, BtJ, ItJ};
	
	velloopjnt1 = VelCOM[ body ] +
		Cross[ AngVel[body], BtJ];
	velloopjnt2 = VelCOM[ inboard ] +
		Cross[ AngVel[inboard], ItJ];
	joint = velloopjnt1 - velloopjnt2;
	angvel = AngVel[body] - AngVel[inboard];
	{vector1,vector2} = If[ PrincipalAxisQ[
		ConvertList[axis,ground] ],	
		ConvertVec[ 
			NullSpace[ {ConvertList[axis,ground]} ], ground],
		ConvertVec[
			NullSpace[ {ConvertList[axis,c]} ] ], ground];	
	constraints = Map[ # == 0 &, Select[ {
		joint . axis,
		joint . vector1,
		joint . vector2,
		angvel . vector1,
		angvel . vector2 }, # =!= 0 & ] ] /. 
		Join[replacesin,replacecos,replace2SC,replaceSC2,
			replaceUs];
		
	Nonholonomic = Join[ Nonholonomic, constraints ];
		
  (* Return the number of constraints introduced *)
  Length[ constraints ]
];


CloseLoop[ body_, inboard_, Fixed, opts___Rule ] := Module[ 
		{BtJ, ItJ, axis, velloopjnt1, velloopjnt2,
		joint, vector1, vector2, angvel, constraints},

	axis = ConvertVec[Axis, body] /. 
		{opts} /. Options[ AddFrame ];
	
 	{BtJ, ItJ} = {BodyToJnt, InbToJnt} /.
		{opts} /. Options[ AddBody ];

	Loop[body] = {inboard, BtJ, ItJ};
	
	velloopjnt1 = VelCOM[ body ] +
		Cross[ AngVel[body], BtJ];
	velloopjnt2 = VelCOM[ inboard ] +
		Cross[ AngVel[inboard], ItJ];
	joint = velloopjnt1 - velloopjnt2;
	angvel = AngVel[body] - AngVel[inboard];
	{vector1,vector2} = If[ PrincipalAxisQ[
		ConvertList[axis,inboard] ],	
		NullSpace[ {ConvertList[axis,inboard]} ],
		NullSpace[ {ConvertList[axis,body]} ] ];
	vector1 = ConvertVec[ vector1, body ];
	vector2 = ConvertVec[ vector2, body ];
	
	constraints = Map[ # == 0 &, Select[ {
		joint . vector1,
		joint . vector2,
		joint . vector3,
		angvel . vector1,
		angvel . vector2,
		angvel . vector3}, # =!= 0 & ] ];
		
	Nonholonomic = Join[ Nonholonomic, constraints ];
		
  (* Return the number of constraints introduced *)
  Length[ constraints ]
];


(* ::Input:: *)
(*End[ ];*)


EndPackage[ ];


(* ::Section:: *)
(*And here's the playground*)


(* ::Input:: *)
(*NewModel[]*)


(* ::Input:: *)
(*AddFrame[a, ground, Hinge, Axis->{0,0,1}, Qdof->1];*)


(* ::Input:: *)
(*AddFrame[a, ground, Hinge, Axis->3]*)


(* ::Input:: *)
(*AddFrame[b, a, Hinge, Axis->{0,0,1}, Qdof->2];*)


(* ::Input:: *)
(*AddFrame[c, ground, Hinge, {0,0,1}];*)


(* ::Input:: *)
(*AddFrame[d, c, Hinge, Axis->{0,0,1}];*)


(* ::Input:: *)
(*AddFrame[e, d, Hinge, Axis->{0,0,1}];*)


(* ::Input:: *)
(*AddFrame[f, b, Hinge, Axis->3];*)


(* ::Input:: *)
(*AddFrame[g, f, Hinge, Axis->3];*)


(* ::Input:: *)
(*?Parents*)


(* ::Input:: *)
(*?Kids*)


(* ::Input:: *)
(*Cross[ a[3] + 2b[3], 3c[2]+3d[1] ]*)


(* ::Section:: *)
(*Transformation Matrix*)


(* ::Special1:: *)
(*Here's how the transformation matrices work:*)
(*Tmtx is the matrix that puts a vector in a frame's Parent coordinates into the frame's coordinates.  For example, for a frame a attached to ground, Tmtx . ground[1] = # a[1] + ... or Transpose[Tmtx] . a[1] = # ground[1] + ...*)
(**)
(*The way we do it is to perform a coordinate transformation that makes the 3 axis point in the direction of the hinge axis.  Then we do an arbitrary coordinate transform about 3, and then transform back into the real world.  *)
(**)
(*What was the hinge becomes 3: if axis is {z1,z2,z3} then {{x1,x2,x3},{y1,y2,y3},{z1,z2,z3}} . axis == {0,0,1}.*)
(**)
(*Then we do the standard 3 axis coordinate transform {{Cos[q],Sin[q],0},{-Sin[q],Cos[q],0},{0,0,1}}*)
(**)
(*Now what was 3 axis becomes hinge: {{x1,y1,z1},{x2,y2,z2},{x3,y3,z3}} . {0,0,1} == {z1,z2,z3}*)
(**)


(* ::Input:: *)
(*GTmtx = Expand[{{x1,y1,z1},{x2,y2,z2},{x3,y3,z3}}.*)
(* {{Cos[q],Sin[q],0},{-Sin[q],Cos[q],0},{0,0,1}}.*)
(* {{x1,x2,x3},{y1,y2,y3},{z1,z2,z3}}]*)


(* ::Special1:: *)
(*We have unknowns in the other 2 axes of our hinge axis world: we only know 3 is the hinge direction, but don't know where 1 and 2 are pointed. Luckily, we know a few things.*)
(*1. It's an orthonormal matrix, so rows and columns dotted against each other give 0.*)
(*2. Cross products: 1 x 2 = 3, 2 x 3 = 1, 3 x 1 = 2.*)
(*3. Norm of rows and columns is 1.*)


(* ::Input:: *)
(*ar1 = AlgebraicRules[ {*)
(*  x1 x2 + y1 y2 + z1 z2 == 0,*)
(*  x1 x3 + y1 y3 + z1 z3 == 0,*)
(*  x2 x3 + y2 y3 + z2 z3 == 0},*)
(*     {x1,x2,x3,y1,y2,y3,z1,z2,z3,Cos[q],Sin[q]} ];*)


(* ::Input:: *)
(*ar2 = AlgebraicRules[ {*)
(*  x2 y3 - x3 y2 == z1,*)
(*  x3 y1 - x1 y3 == z2,*)
(*  x1 y2 - x2 y1 == z3},*)
(*  {x1,x2,x3,y1,y2,y3,z1,z2,z3,Cos[q],Sin[q]} ];*)


(* ::Input:: *)
(* ar3 = AlgebraicRules[ {*)
(*  x1^2 + y1^2 + z1^2 == 1,*)
(*  x2^2 + y2^2 + z2^2 == 1,*)
(*  x3^2 + y3^2 + z3^2 == 1},*)
(*  {x1,x2,x3,y1,y2,y3,z1,z2,z3,Cos[q],Sin[q]} ];*)
(**)


(* ::Input:: *)
(*Collect[ GTmtx /. ar1 /. ar2 /. ar3, {z1,z2,z3} ]*)


(* ::Special1:: *)
(*And by the way, here's the old, uncool way I did this: you could specify the 1, 2, or 3 axis, and it would make a transformation matrix. Problem is, it can't do rotations about arbitrary axes.*)


(* ::Input:: *)
(*AddFrame[frme_,basefrm_,Hinge,coord_?CoordQ,dof_Integer:0 ] :=*)
(*	Module[{newdof = If[ dof != 0, dof, PickDof[dofs] ], mtx},*)
(*	dofs = Union[ Append[dofs, newdof] ];*)
(*	Frames = Union[ Append[Frames, frme] ];*)
(*	*)
(*	(* Now set up the transformation matrix between frames *)*)
(*	mtx = IdentityMatrix[ 3 ];*)
(*	mtx = ReplacePart[mtx, 1, {Mod3[coord],Mod3[coord]}];*)
(*	mtx = ReplacePart[mtx, Cos[q[newdof]],*)
(*  	{{Mod3[coord+1],Mod3[coord+1]},*)
(*  	 {Mod3[coord+2],Mod3[coord+2]}}];*)
(*	mtx = ReplacePart[mtx, -Sin[q[newdof]],*)
(*  	{Mod3[coord+2],Mod3[coord+1]}];*)
(*	mtx = ReplacePart[mtx, Sin[q[newdof]],*)
(*  	{Mod3[coord+1],Mod3[coord+2]}];*)
(*	Tmtx[frme] = mtx;*)
(*	Evaluate[frme][n_?CoordQ] := PV[{1,frme,n}];*)
(*	Parents[frme] = Append[ Parents[basefrm], frme];*)
(*	Kids[frme] = {};*)
(*	Kids[basefrm] = Append[ Kids[basefrm], frme];*)
(*	newdof*)
(*];*)
