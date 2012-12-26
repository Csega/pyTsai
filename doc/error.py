#----------------------------------------------------------
# File error.py
# Simple error dialog
#----------------------------------------------------------
 
import bpy
from bpy.props import *
 
#
#   The error message operator. When invoked, pops up a dialog 
#   window with the given message.
#
class MessageOperator(bpy.types.Operator):
    bl_idname = "error.message"
    bl_label = "Message"
    message = StringProperty()
 
    def execute(self, context):
        self.report({'INFO'}, self.message)
        print(self.message)
        return {'FINISHED'}
 
    def invoke(self, context, event):
        wm = context.window_manager
        return wm.invoke_popup(self, width=400, height=200)
 
    def draw(self, context):
        self.layout.label("ERROR")
        row = self.layout.split(1.0)
        row.prop(self, "message")
 
# Register classes and start scan automatically
bpy.utils.register_class(MessageOperator)

def menuError(str):
    """
	Displays an error string in a pop-up menu.
	"""
    bpy.ops.error.message('INVOKE_DEFAULT', 
        message = str)    
    return

if __name__ == "__main__":
	menuError('valamihiba')