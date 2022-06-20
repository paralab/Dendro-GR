
"""
@author : Milinda Fernando
School of Computing, University of Utah. 
@date: 21 Jan. 2019

@package: Contains utility function to generate C/C++ code, to use with 
Dendro-5.0 octree based PDE solver. 

This is based on [1].
but extended to support templated classes with symbolic Python. 

[1]. https://www.codeproject.com/Articles/571645/Really-simple-Cplusplus-code-generation-in-Python

Copyright 2019 Milinda Fernando

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), 
to deal in the Software without restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense,and/or sell copies
of the Software,and to permit persons to whom the Software is 
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included 
in all copies or substantial portions of the Software. 
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
OTHER DEALINGS IN THE SOFTWARE.
"""
import re
PLACEHOLDER = re.compile('\\$([^\\$]+)\\$')

class Snippet:
	last = None
	def __init__(self, owner, text, postfix):
		self.owner = owner
		if self.owner.last is not None:
			with self.owner.last:
				pass
		self.owner.write("".join(text))
		self.owner.last = self
		self.postfix = postfix
		
	def __enter__(self):
		self.owner.write("{")
		self.owner.current_indent += 1
		self.owner.last = None
		
	def __exit__(self, a, b, c):
		if self.owner.last is not None:
			with self.owner.last:
				pass
		self.owner.current_indent -= 1
		self.owner.write("}" + self.postfix)
		
class Subs:
	def __init__(self, owner, subs):
		self.owner = owner
		self.subs = subs
		
	def __enter__(self):
		self.owner.substack = [self.subs] + self.owner.substack
		
	def __exit__(self, a, b, c):
		self.owner.substack = self.owner.substack[1:]
		

class CodeFile:
	def __init__(self, filename):
		self.current_indent = 0
		self.last = None
		self.out = open(filename,"w")
		self.indent = "\t"
		self.substack = []
		
	def close(self):
		self.out.close()
		self.out = None
	
	def write(self, x, indent=0):
		self.out.write(self.indent * (self.current_indent+indent) + x + "\n")
		
	def format(self, text):
		while True:
			m = PLACEHOLDER.search(text)
			if m is None:
				return text
			s = None
			for sub in self.substack:
				if m.group(1) in sub:
					s = sub[m.group(1)]
					break
			if s is None:
				raise Exception("Substitution '%s' not set." % m.groups(1))
			text = text[:m.start()] + str(s) + text[m.end():]		
		
	def subs(self, **subs):
		return Subs(self, subs)
		
	def __call__(self, text):
		self.write(self.format(text))
		
	def block(self, text, postfix=""):
		return Snippet(self, self.format(text), postfix)

class CppFile(CodeFile):
	def __init__(self, filename):
		CodeFile.__init__(self, filename)
		
	def label(self, text):
		self.write(self.format(text) + ":", -1)
		
__all__ = [ "CppFile", "CodeFile" ]
