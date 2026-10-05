-- On a protocol page (meta `protocol-id`), fill [ ]{.protocol-count} with the number of
-- notebook entries whose front matter lists that id under `protocols:`.
local function count_uses(id)
  local dir = (quarto and quarto.project and quarto.project.directory or ".") .. "/notebook"
  local ok, files = pcall(pandoc.system.list_directory, dir)
  if not ok then return 0 end
  local n = 0
  for _, f in ipairs(files) do
    if f:match("%.q?md$") then
      local fh = io.open(dir .. "/" .. f, "r")
      if fh then
        local head = fh:read(4000) or ""
        fh:close()
        local fm = head:match("^%-%-%-\n(.-)\n%-%-%-") or ""
        local line = fm:match("\nprotocols:%s*([^\n]*)") or fm:match("^protocols:%s*([^\n]*)")
        if line and line:find('"' .. id .. '"', 1, true) then n = n + 1 end
      end
    end
  end
  return n
end

local count = nil

function Meta(m)
  if m["protocol-id"] then count = count_uses(pandoc.utils.stringify(m["protocol-id"])) end
end

function Span(el)
  if count and el.classes:includes("protocol-count") then
    return pandoc.Str(tostring(count))
  end
end

return { { Meta = Meta }, { Span = Span } }
