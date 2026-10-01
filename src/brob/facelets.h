#pragma once
#include <string>
#include "xcube.h"

/*
 *   Kociemba facelet strings for xcube.  For a cube of size S the string has
 *   6 S^2 letters, one per sticker: the faces in the order U R F D L B, each
 *   read row by row as seen in the usual net (U with B at the top, D with F
 *   at the top, the other four with U at the top).  A letter names the face
 *   whose colour the sticker has in the solved cube.  Whitespace is ignored
 *   on input.
 *
 *   xcube_from_facelets fills a cube constructed with the matching size
 *   (Mx = S-2-My, My = (S-2)/2) and returns 1, or returns 0 with a reason in
 *   err if the string is not a cube state: wrong length or letters, a centre
 *   orbit without four of each colour, a missing or repeated corner, middle
 *   edge or wing, centre centres that are not an orientation of the cube,
 *   or a corner twist or middle edge flip no sequence of moves can make.
 */
std::string xcube_to_facelets(const xcube &cube) ;
int xcube_from_facelets(xcube &cube, const std::string &s, std::string &err) ;
