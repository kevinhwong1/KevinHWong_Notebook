-- Adds a "Type · Project · Tools" line under the title of notebook posts and pipelines.
local function str(x) return x and pandoc.utils.stringify(x) or nil end

local function slug(s)
  return (s:lower():gsub("[^%w]+", "-"):gsub("^-+", ""):gsub("-+$", ""))
end

-- Quarto reserves the `project` key, so read it straight from the source file's front matter.
local function project_from_source()
  local f = quarto and quarto.doc and quarto.doc.input_file
  if not f then return nil end
  local fh = io.open(f, "r")
  if not fh then return nil end
  local head = fh:read(3000) or ""
  fh:close()
  local fm = head:match("^%-%-%-\n(.-)\n%-%-%-") or ""
  local p = fm:match('\nproject:%s*"([^"\n]*)"') or fm:match("\nproject:%s*([^\n]*)")
  if p and p ~= "" and p ~= '""' then return p end
  return nil
end

function Pandoc(doc)
  local m = doc.meta
  local parts = {}
  if m.type then
    table.insert(parts, '<span class="pm-type">' .. str(m.type) .. '</span>')
  end
  local p = project_from_source()
  if p then
    table.insert(parts, '<span class="pm-label">Project</span> <a href="../projects.html#' .. slug(p) .. '">' .. p .. '</a>')
  end
  if m.tools then
    local t = {}
    for _, v in ipairs(m.tools) do table.insert(t, '<code>' .. str(v) .. '</code>') end
    if #t > 0 then table.insert(parts, '<span class="pm-label">Tools</span> ' .. table.concat(t, " ")) end
  end
  if m.protocols then
    local root = (quarto and quarto.project and quarto.project.directory) or "."
    local links = {}
    for _, v in ipairs(m.protocols) do
      local id = str(v)
      local title = id
      local fh = io.open(root .. "/protocols/" .. id .. ".qmd", "r")
      if fh then
        local head = fh:read(2000) or ""
        fh:close()
        title = head:match('\ntitle:%s*"([^"\n]+)"') or title
      end
      table.insert(links, '<a href="../protocols/' .. id .. '.html">' .. title .. '</a>')
    end
    if #links > 0 then
      table.insert(parts, '<span class="pm-label">Protocol</span> ' .. table.concat(links, ", "))
    end
  end
  if m.pipeline then
    table.insert(parts, '<span class="pm-label">Pipeline</span> <a href="' .. str(m.pipeline) .. '">view</a>')
  end
  if #parts > 0 then
    local html = '<div class="post-meta">' .. table.concat(parts, '<span class="pm-sep">·</span>') .. '</div>'
    table.insert(doc.blocks, 1, pandoc.RawBlock("html", html))
  end
  return doc
end
