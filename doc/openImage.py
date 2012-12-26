import bpy
 
class OpenImage(bpy.types.Operator):
    bl_idname = "dialog.open_image"  # this is important since its how bpy.ops.dialog.open_file is constructed
    bl_label = "Open Image"
    filepath = bpy.props.StringProperty(subtype="FILE_PATH")
    x = bpy.props.IntProperty()
    y = bpy.props.IntProperty()
    
    def invoke(self, context, event):
        wm = context.window_manager.fileselect_add(self)
        self.x = event.mouse_x
        self.y = event.mouse_y
        return {'RUNNING_MODAL'}
    
    def execute(self, context):
        # rather then printing, use the report function,
        # this way the messag appiers in the header,
        #self.report({'INFO'}, "Mouse coords are %d %d" % (self.x, self.y))
        #self.report({'INFO'},"The file path is %s" % (self.filepath))
        bpy.ops.image.open(filepath=self.filepath)
		bpy.data.texts.new('image_filepath')
		bpy.data.texts['image_filepath'].write(self.filepath)
        #self.report({'INFO'},"The file path is %s" % (self.filepath))
        #self.report({'INFO'},"The images are %s" % (bpy.data.images.keys()))
        return {'FINISHED'}
    
bpy.utils.register_class(OpenImage)
#bpy.ops.dialog.open_image('INVOKE_DEFAULT')
