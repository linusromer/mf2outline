#!/usr/bin/env python
import os,sys,fontforge,glob,subprocess,tempfile,shutil
if __name__ == "__main__":			
	mffile = os.path.abspath(sys.argv[1])
	tempdir = tempfile.mkdtemp()
	font = fontforge.font()
	subprocess.call(
	['mpost',
	'&mfplain',
	'\mode=localfont;',
	'mag:=100.375;',
	'outputtemplate:="%c.eps";',
	'input %s;' % mffile,
	'bye'],
	stdout=subprocess.PIPE, 
	stderr=subprocess.PIPE,
	cwd=tempdir
	)
	glyph_files = glob.glob(os.path.join(tempdir, "*.eps"))
	for eps in glyph_files:
		code = int(os.path.splitext(os.path.basename(eps))[0])
		glyph = font.createChar(code)
		glyph.importOutlines(eps, ("toobigwarn", "correctdir"))
	font.generate("font.otf")
	shutil.rmtree(tempdir)
	exit(0)
